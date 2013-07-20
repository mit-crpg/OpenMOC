#include "Cell.h"


short int Cell::_n = 0;

static int auto_id = 10000;


/**
 * @brief Returns an auto-generated unique cell ID.
 * @details This method is intended as a utility mehtod for user's writing
 *          OpenMOC input files. The method makes use of a static cell
 *          ID which is incremented each time the method is called to enable
 *          unique generation of monotonically increasing IDs. The method's
 *          first ID begins at 10000. Hence, user-defined cell IDs greater
 *          than or equal to 10000 is prohibited.
 */
int cell_id() {
    int id = auto_id;
    auto_id++;
    return id;
}


/**
 * @brief Default constructor used in rings/sectors subdivision of cells.
 */
Cell::Cell() { }


/**
 * @brief Constructor sets the unique and user-specifed IDs for this cell.
 * @param id the user-specified cell ID
 * @param universe the ID of the universe within which this cell resides
 */
Cell::Cell(short int universe, short int id) {

    /* If the user did not define an optional ID, create one */
    if (id == 0)
        _id = surf_id();
    else if (id >= surf_id())
        log_printf(ERROR, "Unable to set the ID of a cell to %d since "
		 "cell IDs greater than or equal to 10000 is probibited "
		 "by OpenMOC.", id);
    /* Use the user-defined ID */
    else
        _id = id;

    _uid = _n;
    _n++;
    _universe = universe;       
}


/**
 * @brief Destructor frees all surfaces making up cell.
 */
Cell::~Cell() {
    _surfaces.clear();
}


/**
 * @brief Return the cell's unique ID.
 * @return the cell's unique ID
 */
short int Cell::getUid() const {
    return _uid;
}


/**
 * @brief Return the cell's user-specified ID.
 * @return the cell's user-specified ID
 */
short int Cell::getId() const {
    return _id;
}


/**
 * @brief Return the cell type (FILL or MATERIAL).
 * @return the cell type
 */
cellType Cell::getType() const {
    return _cell_type;
}


/**
 * @brief Return the ID of the universe within which this cell resides.
 * @return the universe ID
 */
short int Cell::getUniverse() const {
    return _universe;
}


/**
 * @brief Return the number of surfaces in the cell.
 * @return the number of surfaces
 */
short int Cell::getNumSurfaces() const {
    return this->_surfaces.size();
}


/**
 * @brief Return the hashtable of surfaces IDs and surface pointers for all 
 *        surfaces making up the cell.
 * @return map of surface IDs and surface pointers
 */
std::map<short int,Surface*> Cell::getSurfaces() const {
    return _surfaces;
}


/**
 * @brief Set the ID for teh universe within whic this cell resides.
 * @param universe the universe's user-specified ID
 */
void Cell::setUniverse(short int universe) {
    _universe = universe;
}


/**
 * @brief Insert the a surface into this cell's container.
 * @param halfspace the surface halfspace (positive/negative for surface side)
 * @param surface a pointer to the surface
 */
void Cell::addSurface(short int halfspace, Surface* surface) {
    if (halfspace != -1 && halfspace != +1)
        log_printf(ERROR, "Unable to add surface %d to cell %d since the "
                    "halfspace %d is not -1 or 1", surface->getId(), 
                                                    _id, halfspace);

    _surfaces.insert(std::pair<short int, Surface*>(halfspace*surface->getId(),
                                                                     surface));
}


/**
 * @brief Registers a surface pointer with the Cell's surfaces map. 
 * @details This method is used by the geometry whenever a new cell is 
 *          added to it.
 * @param surface the surface pointer
 */
void Cell::setSurfacePointer(Surface* surface) {

    /* if _surfaces does not contain this surface id throw an error */
    if (_surfaces.find(surface->getId()) == _surfaces.end() &&
          _surfaces.find(-surface->getId()) == _surfaces.end())

      log_printf(WARNING, "Unable to set surface pointer for cell id = %d "
                 "for surface id = %d since cell does not contain this surface",
                 _id, surface->getId());

    try{
        /* If the cell contains the positive side of the surface */
        if (_surfaces.find(surface->getId()) != _surfaces.end())
            _surfaces[surface->getId()] = surface;

        /* If the cell contains the negative side of the surface */
        else
            _surfaces[-1*surface->getId()] = surface;

        log_printf(INFO, "Set the surface pointer for cell id = %d for "
                   "surface id = %d", _id, surface->getId());
    }
    catch (std::exception &e) {
        log_printf(ERROR, "Unable to add surface with id = %d to cell with "
		 "id = %d. Backtrace:\n%s", surface, _id, e.what());
    }
}


/**
 * @brief Determines whether a point is contained inside a cell. 
 * @details Queries each surface inside the cell to determine if the point
 *          particle is on the same side of the surface. This point is 
 *          only inside the cell if it is on the same side of every surface 
 *          in the cell.
 * @param point a pointer to a point
 */
bool Cell::cellContainsPoint(Point* point) {

    /* Loop over all surfaces inside the cell */
    std::map<short int, Surface*>::iterator iter;
    for (iter = _surfaces.begin(); iter != _surfaces.end(); ++iter) {
      
        /* If the surface evaluated point is not the same sign as the surface
         * or within a threshold, return false */
        if (iter->second->evaluate(point) * iter->first < -ON_SURFACE_THRESH)
            return false;
    }

    return true;
}


/**
 * @brief Determines whether a point is contained inside a cell. 
 * @details Queries each surface inside the cell to determine if the particle 
 *          is on the same side of the surface. This particle is only inside 
 *          the cell if it is on the same side of every surface in the cell.
 * @param coords a pointer to a localcoord
 */
bool Cell::cellContainsCoords(LocalCoords* coords) {
    return this->cellContainsPoint(coords->getPoint());
}


/**
 * @brief Computes the minimum distance to a surface from a point with a given
 *        trajectory at a certain angle. 
 * @details If the trajectory will not intersect any of the surfaces in the cell*           returns INFINITY.
 * @param point the point of interest
 * @param angle the angle of the trajectory (in radians from 0 to 2*PI)
 * @param min_intersection a pointer to the intersection point that is found
 */
double Cell::minSurfaceDist(Point* point, double angle,
                            Point* min_intersection) {

    double min_dist = INFINITY; 
    double d;
    Point intersection;

    std::map<short int, Surface*>::iterator iter;
    
    /* Loop over all of the cell's surfaces */
    for (iter = _surfaces.begin(); iter != _surfaces.end(); ++iter) {

        /* Find the minimum distance from this surface to this point */
        d = iter->second->getMinDistance(point, angle, &intersection);

        /* If the distance to cell is less than current min distance, update */
        if (d < min_dist) {
            min_dist = d;
            min_intersection->setX(intersection.getX());
            min_intersection->setY(intersection.getY());
        }
    }

    return min_dist;
}


/**
 * Constructor sets the user-specified and unique IDs for this cellbasic.
 * @param universe the ID for the universe within which this cellbasic resides
 * @param material the ID for the material which fills this cellbasic
 * @param rings the number of equal volume rings to divide this cell into 
 *        (the default is zero)
 * @param sectors the number of angular sectors to divide this cell into
 *        (the default is zero)
 * @param id the user-specified cell ID
 */
CellBasic::CellBasic(short int universe, short int material, int rings, 
		     int sectors, short int id): Cell(universe, id) {
    _cell_type = MATERIAL;
    _material = material;
    setNumRings(rings);
    setNumSectors(sectors);
}



/**
 * @brief Return the ID of the material filling the cellbasic.
 * @return the material's ID
 */
short int CellBasic::getMaterial() const {
    return _material;
}


/**
 * @brief Return the number of rings in the cell.
 * @return the number of rings
 */
short int CellBasic::getNumRings() {
    return this->_num_rings;
}


/**
 * @brief Return the number of sectors in the cell.
 * @return the number of sectors
 */
short int CellBasic::getNumSectors() {
    return this->_num_sectors;
}


/**
 * @brief Return the number of flat source regions in this cellbasic. 
 * @details This method is used when the geometry recursively constructs flat
 *           source regions. By definition, cellbasic's are the lowest level
 *           entity in the geometry and thus only have one flat source region
 *           within them, so this method always returns 1.
 * @return the number of FSRs in this cell
 */
int CellBasic::getNumFSRs() {
    return 1;
}


/**
 * @brief Set the cell's number of rings.
 * @param num_rings the number of rings in this cell
 */
void CellBasic::setNumRings(short int num_rings) {
    if (num_rings < 0)
        log_printf(ERROR, "Unable to give %d rings to cell %d since this is "
		 "a negative number", num_rings, _id);

    _num_rings = num_rings;
}


/**
 * @brief Set the cell's number of sectors.
 * @param num_sectors the number of sectors in this cell
 */
void CellBasic::setNumSectors(short int num_sectors) {
    if (num_sectors < 0)
        log_printf(ERROR, "Unable to give %d sectors to cell %d since this is "
		 "a negative number", num_sectors, _id);

    if (num_sectors == 1)
        _num_sectors = 0;
    else
        _num_sectors = num_sectors;
}


/**
 * @brief Create a duplicate of the cellbasic.
 * @return a pointer to the clone
 */
CellBasic* CellBasic::clone() {

    /* Construct new cell */
    CellBasic* new_cell = new CellBasic(_universe, _material, 
					_num_rings, _num_sectors);

    /* Loop over all of this cell's surfaces and add them to the clone */
    std::map<short int, Surface*>::iterator iter;
    int halfspace;
    Surface* surface;

    for (iter = _surfaces.begin(); iter != _surfaces.end(); ++iter) {
        halfspace = iter->first / abs(iter->first);
        surface = iter->second;
        new_cell->addSurface(halfspace, surface);
    }

    return new_cell;
}


void CellBasic::sectorize() { 

    if (_num_sectors == 0)
        return;

    /* Figure out the angle for each sector */
    double* azim_angles = new double[_num_sectors];
    double delta_azim = 2. * M_PI / _num_sectors;
    double A, B;

    std::vector<Plane*> planes;
    std::vector<Plane*>::iterator iter1;

    log_printf(DEBUG, "Sectorizing cell %d with %d sectors",_id, _num_sectors);

    for (int i=0; i < _num_sectors; i++) {
        azim_angles[i] = i * delta_azim;
	A = cos(azim_angles[i]);
	B = sin(azim_angles[i]);
	Plane* plane = new Plane(A, B, 0.);
	planes.push_back(plane);
	log_printf(DEBUG, "Created sector plane id = %d, angle = %f, A = %f, "
		   "B = %f", i, azim_angles[i], A, B);
    }

    /* Create sectors using disjoint halfspaces of pairing planes */
    for (int i=0; i < _num_sectors; i++) {

        /* Create new CellBasic clones */
        CellBasic* sector = clone();

        sector->setNumSectors(0);
        sector->setNumRings(0);

        log_printf(DEBUG, "Creating a new sector cell with %d for cell %d", 
		   sector->getId(), _id);

        /* Add new bounding planar surfaces to the clone */
        sector->addSurface(+1, planes.at(i));

	if (_num_sectors != 2) {
	    if (i+1 < _num_sectors)
	        sector->addSurface(-1, planes.at(i+1));
	    else
	        sector->addSurface(-1, planes.at(0));
	}

	_sectors.push_back(sector);
    }

    _subcells.clear();
    _subcells.insert(_subcells.end(), _sectors.begin(), _sectors.end());

    delete [] azim_angles;
}


void CellBasic::ringify() {

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
    short int halfspace1 = 0;
    short int halfspace2 = 0;
    std::vector<Circle*> circles;

    /* See if the cell contains 1 or 2 surfaces */
    std::map<short int, Surface*>::iterator iter1;
    for (iter1=_surfaces.begin(); iter1 != _surfaces.end(); ++iter1) {
        if ((*iter1).second->getSurfaceType() == CIRCLE) {
	    short int halfspace = (*iter1).first / (*iter1).second->getId();
	    Circle* circle = static_cast<Circle*>((*iter1).second);

	    /* Outermost bounding circle */
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
        log_printf(ERROR, "Unable to ringify cell %d since it does not "
		   "contain any CIRCLE type surface(s)", _id);
    if (num_circles > 2)
        log_printf(NORMAL, "Unable to ringify cell %d since it "
		   "contains more than 2 CIRCLE surfaces", _id);
    if (x1 != x2 && num_circles == 2)
        log_printf(ERROR, "Unable to ringify cell %d since it contains "
		 "circle %d centered at x=%f and circle %d at x=%f. "
		   "Both circles must have the same center.", 
		   _id, circle1->getId(), x1, circle2->getId(), x2);
    if (y1 != y2 && num_circles == 2)
        log_printf(ERROR, "Unable to ringify cell %d since it contains "
		 "circle %d centered at y=%f and circle %d at y=%f. "
		   "Both circles must have the same center.", 
		   _id, circle1->getId(), y1, circle2->getId(), y2);
    if (circle1 == NULL && circle2 != NULL)
        log_printf(ERROR, "Unable to ringify cell %d since it only contains "
		   "the positive halfpsace of circle %d. Rings can only be "
		   "created for cells on the interior (negative halfspace) "
		   "of a circle surface.", _id, circle2->getId());
    if (radius1 <= radius2)
        log_printf(ERROR, "Unable to ringify cell %d since it contains 2 "
		 "disjoint CIRCLE surfaces: halspace %d for circle %d "
		   "and halfspace %d for circle %d. Switch the signs of "
		   "the 2 halfspaces for each surface.", _id, halfspace1,
		   circle1->getId(), halfspace2, circle2->getId());

    /* Compute the area to fill with each ring */
    double area = M_PI * fabs(radius1*radius1 - radius2*radius2) / _num_rings;

    /* Generate successively smaller circle surfaces */
    for (int i=0; i < _num_rings-1; i++) {
        radius2 = sqrt(radius1*radius1 - (area / M_PI));
	Circle* circle = new Circle(x1, y1, radius1);
	circles.push_back(circle);
	radius1 = radius2;
    }
 
    /* Store smallest, innermost circle */
    Circle* circle = new Circle(x1, y1, radius1);
    circles.push_back(circle);

    /* Loop over circles and */
    std::vector<Circle*>::iterator iter2;
    std::vector<CellBasic*>::iterator iter3;

    for (iter2 = circles.begin(); iter2 != circles.end(); ++iter2) {

        /* Create circles for each of the sectorized cells */
	if (_sectors.size() != 0) {
	    for (iter3 = _sectors.begin(); iter3 != _sectors.end(); ++iter3) {
	        log_printf(DEBUG, "Creating a new ring in sector cell %d",
			 (*iter3)->getId());

                /* Create a new CellBasic clone */
	        CellBasic* ring = (*iter3)->clone();
	        ring->setNumSectors(0);
	        ring->setNumRings(0);

	        /* Add new bounding circle surfaces to the clone */
	        ring->addSurface(-1, (*iter2));

	        /* Look ahead and check if we have an inner circle to add */
	        if (iter2+1 == circles.end()) {
	            _rings.push_back(ring);
		    continue;
	        }
	        else
	            ring->addSurface(+1, *(iter2+1));

		_rings.push_back(ring);
	    }
	}

	/* Create circles for this un-sectorized cell */
	else {
	    log_printf(DEBUG, "Creating new ring in un-sectorized cell %d",_id);

            /* Create a new CellBasic clone */
	    CellBasic* ring = clone();
	    ring->setNumSectors(0);
	    ring->setNumRings(0);

	    /* Add new bounding circle surfaces to the clone */
	    ring->addSurface(-1, (*iter2));
	    
	    /* Look ahead and check if we have an inner circle to add */
	    if (iter2+1 == circles.end()) {
	        _rings.push_back(ring);
		break;
	    }
	    else
	        ring->addSurface(+1, *(iter2+1));

	    _rings.push_back(ring);
	}
    }

    _subcells.clear();
    _subcells.insert(_subcells.end(), _rings.begin(), _rings.end());
}


/**
 * @brief Subdivides a cells into rings and sectors.
 * @details This method uses the Cell's clone method to produce
 *          a vector of new cells, each representing a subdivision
 *          of this cell into rings and sectors.
 * @return a vector of Cell pointers to the new subdivided cells
 */
std::vector<CellBasic*> CellBasic::subdivideCell() {
    sectorize();
    ringify();
    return _subcells;
}

/**
 * @brief Convert this cellbasic's attributes to a string format.
 * @return a character array of this cellbasic's attributes
 */
std::string CellBasic::toString() {

    std::stringstream string;

    string << "Cell id = " << _id 
           << ", type = MATERIAL, material id = " << _material 
           << ", universe = " << _universe  
           << ", num_surfaces = " << getNumSurfaces() 
           << ", num of rings = " << _num_rings 
           << ", num of sectors = " << _num_sectors;
    
    /* Append each of the surface ids to the string */
    std::map<short int, Surface*>::iterator iter;
    string << ", surface ids = ";
    for (iter = _surfaces.begin(); iter != _surfaces.end(); ++iter)
        string << iter->first << ", ";

    return string.str();
}


/**
 * @brief Prints a string representation of all of the cellbasic's objects to
 *        the console.
 */
void CellBasic::printString() {
    log_printf(RESULT, toString().c_str());
}


/**
 *  @brief CellFill constructor
 *  @param id the user-specified cell ID
 *  @param universe the ID of the universe within which this cell resides
 *  @param universe_fill the ID of the universe filling this cell
 */
CellFill::CellFill(short int universe, short int universe_fill, short int id):
  Cell(universe, id) {
        _cell_type = FILL;
        _universe_fill.first = universe_fill;
}


/**
 * @brief Return the ID of the universe filling this cell.
 * @return the universe's id
 */
short int CellFill::getUniverseFillId() const {
    return _universe_fill.first;
}


/**
 * @brief Return a pointer to the universe filling this cell.
 * @return the universe pointer
 */
Universe* CellFill::getUniverseFill() const {
    return _universe_fill.second;
}


/**
 * @brief Return the number of flat source regions in this cellfill. 
 * @details This method is used when the geometry recursively constructs flat
 *          source regions.
 * @return the number of FSRs in this cellfill
 */
int CellFill::getNumFSRs() {
    Universe *univ = getUniverseFill();
    if (univ->getType() == SIMPLE)
        return univ->computeFSRMaps();
    else
        return static_cast<Lattice*>(univ)->computeFSRMaps();
}


/**
 * @brief Set the ID of the universe filling this cellfill.
 * @param universe_fill the universe's ID
 */
void CellFill::setUniverseFill(short int universe_fill) {
    _universe_fill.first = universe_fill;
}


/**
 * @brief Set a pointer to the universe filling this cellfill.
 * @param universe the universe's pointer
 */
void CellFill::setUniverseFillPointer(Universe* universe) {
    _universe_fill.second = universe;
}


/**
 * @brief Convert this cellfill's attributes to a string format.
 * @return a character array of this cell's attributes
 */
std::string CellFill::toString() {

    std::stringstream string;

    string << "Cell id = " << _id << ", type = FILL, universe_fill = " <<
        _universe_fill.first << ", universe = " << _universe <<
      ", num_surfaces = " << getNumSurfaces();

    /** Add the IDs for the surfaces in this cell */
    std::map<short int, Surface*>::iterator iter;
    string << ", surface ids = ";
    for (iter = _surfaces.begin(); iter != _surfaces.end(); ++iter)
        string << iter->first << ", ";

    return string.str();
}


/**
 * @brief Prints a string representation of all of the cellfill's objects to
 *        the console.
 */
void CellFill::printString() {
    log_printf(RESULT, toString().c_str());
}
