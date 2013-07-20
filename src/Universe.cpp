#include "Universe.h"


short int Universe::_n = 0;


/**
 * @brief Constructor assigns a unique and user-specified ID for the universe.
 * @param id the user-specified universe id
 */
Universe::Universe(const short int id) {
    _uid = _n;
    _id = id;
    _n++;
    _type = SIMPLE;
}


/**
 * @brief Destructor clears the cell pointers container.
 */
Universe::~Universe() {
    _cells.clear();
}


/**
 * @brief Returns the universe's unique ID.
 * @return the universe's unique ID.
 */
short int Universe::getUid() const {
    return _uid;
}

/**
 * Return the user-specified ID for this universe.
 * @return the user-specified universe ID
 */
short int Universe::getId() const {
    return _id;
}


/**
 * @brief Return the universe type (SIMPLE or LATTICE).
 * @return the universe type
 */
universeType Universe::getType() {
    return _type;
}

/**
 * @brief Return the number of cells in this universe.
 * @return the number of cells
 */
short int Universe::getNumCells() const {
    return _cells.size();
}


/**
 * @brief Return a pointer to the origin for this cell (in global coordinates).
 * @return the origin of the cell
 */
Point* Universe::getOrigin() {
    return &_origin;
}


/**
 * @brief Returns the local ID for the FSR representing a cell in this universe.
 * @details This method is used when constructing an ID for a FSR.
 * @param cell_id the ID of the cell of interest
 */
int Universe::getFSR(short int cell_id) {

    if (_cells.find(cell_id) == _cells.end())
        log_printf(ERROR, "Tried to find FSR id for cell with id = %d in "
		   " universe with id = %d but no cell exists", cell_id, _id);

    return _region_map.at(cell_id);
}


/**
 * @brief Return the vector of cell pointers in this universe.
 * @return vector of cell ids
 */
std::map<short int, Cell*> Universe::getCells() const {
    return _cells;
}


/**
 * @brief Adds a cell to this universe.
 * @details Stores the user-specified cell ID and cell pointer in a hash
 *          table along with all of other cells added to this universe.
 * @param cell the cell id
 */
void Universe::addCell(Cell* cell) {

    try {
        _cells.insert(std::pair<short int, Cell*>(cell->getId(), cell));
	log_printf(INFO, "Added cell with id = %d to universe with id = %d",
		   cell->getId(), _id);
	}
    catch (std::exception &e) {
        log_printf(ERROR, "Unable to add cell with id = %d to universe with"
		   " id = %d. Backtrace:\n%s", cell, _id, e.what());
    }
}


/**
 * @brief Sets the universe type to SIMPLE or LATTICE.
 * @param type the universe type
 */
void Universe::setType(universeType type) {
    _type = type;
}


/**
 * @brief Set the origin in global coordinates for this universe.
 * @param origin a pointer to the origin
 */
void Universe::setOrigin(Point* origin) {
    _origin.setX(origin->getX());
    _origin.setY(origin->getY());
}


/**
 * @brief Finds the cell for which a localcoords object resides.
 * @details Finds the cell that a localcoords object is located inside by 
 *          checking each of this universe's cells. Returns NULL if the 
 *          localcoords is not in any of the cells.
 * @param coords a pointer to the localcoords of interest
 * @param universes a container of all of the universes passed in by geometry
 * @return a pointer the cell where the localcoords is located
 */
Cell* Universe::findCell(LocalCoords* coords,
			 std::map<short int, Universe*> universes) {

    Cell* return_cell = NULL;
    std::map<short int, Cell*>::iterator iter;

    /* Sets the localcoord type to UNIV at this level */
    coords->setType(UNIV);

    /* Loop over all cells in this universe */
    for (iter = _cells.begin(); iter != _cells.end(); ++iter) {
        Cell* cell = iter->second;

	if (cell->cellContainsCoords(coords)) {
	  
	    /* Set the cell on this level */
	    coords->setCell(cell->getId());

	    /* MATERIAL type cell - lowest level, terminate search for cell */
	    if (cell->getType() == MATERIAL) {
	        coords->setCell(cell->getId());
		return_cell = cell;
		return return_cell;
	    }

	    /* FILL type cell - cell contains a universe at a lower level
	     * Update coords to next level and continue search */
	    else if (cell->getType() == FILL) {

	        LocalCoords* next_coords;

		if (coords->getNext() == NULL)
		    next_coords = new LocalCoords(coords->getX(),
						  coords->getY());
		else
		    next_coords = coords->getNext();
		
		CellFill* cell_fill = static_cast<CellFill*>(cell);
		short int universe_id = cell_fill->getUniverseFillId();
		next_coords->setUniverse(universe_id);
		Universe* univ = universes.at(universe_id);
		coords->setCell(cell->getId());

		coords->setNext(next_coords);
		next_coords->setPrev(coords);
                if (univ->getType() == SIMPLE)
  		    return univ->findCell(next_coords, universes);
                else
		    return static_cast<Lattice*>(univ)
		      ->findCell(next_coords, universes);
	    }
	}
    }

    return return_cell;
}


/**
 * @brief Convert the member attributes of this universe to a character array
 * @return a character array reprsenting the universe
 */
std::string Universe::toString() {
    std::stringstream string;
    std::map<short int, Cell*>::iterator iter;

    string << "Universe id = " << _id << ", type = ";
    if (_type == SIMPLE)
        string << "SIMPLE";
    else
        string << "LATTICE";

    string << ", num cells = " << _cells.size() << ", cell ids = ";
    
    for (iter = _cells.begin(); iter != _cells.end(); ++iter)
        string << iter->first << ", ";

    return string.str();
}


/**
 * @brief Prints a string representation of all of the universe's objects to
 *        the console.
 */
void Universe::printString() {
    log_printf(RESULT, toString().c_str());
}


/**
 * @brief Compute the FSR maps for this universe and return the number of
 *        FSRs inside the universe.
 * @details The FSR map is simply a hash table of the local FSR IDs for each
 *          of the cells making up this universe. This is used when 
 *          constructing a global FSR ID for a lattice which may contain
 *          several repeating versions of this universe.
 * @return the number of FSRs in the universe
 */
int Universe::computeFSRMaps() {

    /* initialize a counter count */
    std::map<short int, Cell*>::iterator iter;
    int count = 0;
    
    /* loop over cells in the universe to set the map and update count */
    for (iter = _cells.begin(); iter != _cells.end(); ++iter) {
        _region_map.insert(std::pair<int, int>(iter->first, count));
	count += iter->second->getNumFSRs();
    }

    return count;
}


/**
 * @brief Subdivides all of the cells within this universe into rings
 *        and angular sectors.
 */
void Universe::subdivideCells() {

    log_printf(DEBUG, "Subdividing cells for universe %d", _id);

    std::map<short int, Cell*>::iterator iter1;

    while (iter1 != _cells.end()) {

        for (iter1 = _cells.begin(); iter1 != _cells.end(); ++iter1) {

	    if ((*iter1).second->getType() == MATERIAL) {
	        CellBasic* cell = static_cast<CellBasic*>((*iter1).second);

		if (cell->getNumRings() > 0 || cell->getNumSectors() > 0) {

		    std::vector<CellBasic*> newcells = cell->subdivideCell();
		    
		    log_printf(DEBUG, "cell %d in universe %d has %d subcells",
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
 * @brief Constructor sets the user-specified and unique IDs for this lattice.
 * @param id the user-specified lattice (universe) ID
 * @param width_x the width of the lattice cells along x
 * @param width_y the width of the lattice cells along y
 */
Lattice::Lattice(const short int id, double width_x, double width_y): 
    Universe(id) {

    _width_x = width_x;
    _width_y = width_y;
    _type = LATTICE;

    /* Default number of lattice cells along each dimension */
    _num_y = 0;
    _num_x = 0;
}


/**
 * @brief Destructor clears memory for all of universes pointers.
 */
Lattice::~Lattice() {
  
  for (int i=0; i < _num_x; i++){
      _universes.at(i).clear();
  }

  _universes.clear();
}


/**
 * @brief Return the number of lattice cells along the x-axis.
 * @return the number of lattice cells along x
 */
short int Lattice::getNumX() const {
    return _num_x;
}


/**
 * @brief Return the number of lattice cells along the y-axis
 * @return the number of lattice cells along y
 */
short int Lattice::getNumY() const {
    return _num_y;
}


/**
 * @brief Return the origin of the lattice.
 * @return a pointer to the origin of the lattice
 */
Point* Lattice::getOrigin() {
    return &_origin;
}


/**
 * @brief Return a 2D vector array of the universes in the lattice.
 * @return 2D vector of universes
 */
std::vector< std::vector< std::pair<short int, Universe*> > > 
Lattice::getUniverses() const {
    return _universes;
}


/**
 * @brief Returns a pointer to the universe within a specific lattice cell.
 * @param lattice_x the x index to the lattice cell
 * @param lattice_y the y index to the lattice cell
 * @return pointer to a universe filling the lattice cell
 */
Universe* Lattice::getUniverse(short int lattice_x,short int lattice_y) const {

    /* Checks that lattice indices are within the bounds of the lattice */
    if (lattice_x > _num_x || lattice_y > _num_y)
        log_printf(ERROR, "Cannot retrieve universe from lattice id = %d: "
		   "Index out of bounds: Tried to access cell x = %d, y = %d "
		   "but bounds are x = %d, y = %d", _id, lattice_x, lattice_y,
		   _num_x, _num_y);

    return _universes.at(lattice_y).at(lattice_x).second;
}


/**
 * @brief Return the width of the lattice along the x-axis.
 * @return the width of the lattice cells along x
 */
double Lattice::getWidthX() const {
    return _width_x;
}


/**
 * @brief Return the width of the lattice along the y-axis.
 * @return the width of the lattice cells along y
 */
double Lattice::getWidthY() const {
    return _width_y;
}


/**
 * @brief Return the id of a flat source region base index (smallest FSR 
 *        region id within a specific lattice cell)
 * @param lat_x the x index of the lattice cell
 * @param lat_y the y index of the lattice cell
 * @return the index of the base FSR
 */
int Lattice::getFSR(short int lat_x, short int lat_y) {

    /* Check if lattice indices are out of bounds */
    if (lat_x > _num_x || lat_y > _num_y)
        log_printf(ERROR, "Tried to access FSR map of lattice id = %d, but "
		   "indices lat_x = %d and lat_y = %d were out of bounds", _id,
		   lat_x, lat_y);

    return _region_map[lat_y][lat_x].second;
}


/**
 * @brief Sets the pointer to a universe filling one of this lattice's
 *        lattice cells.
 * @param universe pointer to the universe 
 */
void Lattice::setUniversePointer(Universe* universe) {

    /* Check that _surfaces contains this surface id and delete the id
     *  otherwise throw an error */
    short int universe_id = universe->getId();
    bool universe_not_found = true;

    /* Loop over all universes in lattice */
    for (int i = 0; i < _num_y; i++) {
        for (int j = 0; j< _num_x; j++) {
	    /* Checks if the universe id matches what the lattice's
	     * universes container expects */
	    if (_universes.at(i).at(j).first == universe_id) {
	        _universes[i][j].second = universe;
		universe_not_found = false;
	    }
	}
    }

    if (universe_not_found)
        log_printf(WARNING, "Tried to set the universe pointer for "
		   "lattice id = %d for universe id = %d but the lattice "
		   "does not contain the universe", _id, universe_id);
    else
        log_printf(INFO, "Set the universe pointer for lattice "
		   "id = %d for universe id = %d", _id, universe_id);

    return;
}


/**
 * @brief Sets the arrary of universe IDs filling each lattice cell.
 * @param universes the array of universe Ids for each lattice cell
 * @param num_x the number of lattice cells along x
 * @param num_y the number of lattice cells along y
 */
void Lattice::setLatticeCells(int num_x, int num_y, short* universes) {

    _num_x = num_x;
    _num_y = num_y;

    /* Set origin to lower left corner with (x=0, y=0) in center of lattice */
    _origin.setX(-_width_x*_num_x/2.0);
    _origin.setY(-_width_y*_num_y/2.0);

    Universe* empty_universe_pointer = NULL;
    
    /* The parser gives the lattice cells in row major order starting from the
     * upper left corner. This double loop reorders the lattice cells from the
     * to start from the lower left corner */
    for (int i = 0; i < _num_y; i++) {
        _universes.push_back(std::vector< std::pair<short int, Universe*> >());
	    for (int j = 0; j< _num_x; j++){
	        _universes.at(i).push_back(std::pair<short int, Universe*>
			           (universes[(_num_y-1-i)*_num_x + j], 
                                                    empty_universe_pointer));
    	}
    }

    /* intialize _region_map */
    for (int i = 0; i < _num_y; i++) {
        _region_map.push_back(std::vector< std::pair<int, int> >());
   	    for (int j = 0; j < _num_x; j++) {
	        _region_map.at(i).push_back(std::pair<int, int>());
    	}
    }
}


/**
 * @brief Checks if a point is within the bounds of a lattice.
 * @param point a pointer to the point of interest
 * @return true if the point is in the bounds, false if not
 */
bool Lattice::withinBounds(Point* point) {

    /* Computes the lattice bounds */
    double bound_x = _num_x/2.0 * _width_x;
    double bound_y = _num_y/2.0 * _width_y;

    double x = point->getX();
    double y = point->getY();

    /* If the point is outside the x bounds */
    if (x > bound_x || x < -1*bound_x)
        return false;
    /* If the point is outside the y boounds */
    else if (y > bound_y || y < -1*bound_y)
        return false;
    /* If the point is within the bounds */
    else
        return true;
}


/**
 * @brief Finds the cell within this lattice that a localcoords is in. 
 * @details This method first find the lattice cell, then searches the 
 *          universe inside that lattice cell. If localcoords is outside 
 *          the bounds of the lattice, this method will return NULL.
 * @param coords the localcoords of interest
 * @param universes a map of all universes passed in from the geometry
 * @return a pointer to the cell this localcoord is in or NULL
 */
Cell* Lattice::findCell(LocalCoords* coords,
			std::map<short int, Universe*> universes) {

    /* Set the localcoord to be a LAT type at this level */
    coords->setType(LAT);

    /* Compute the x and y indices for the lattice cell this coord is in */
    short int lat_x = (int)floor((coords->getX() - _origin.getX()) / _width_x);
    short int lat_y = (int)floor((coords->getY() - _origin.getY()) / _width_y);

    /* Check if the localcoord is on the lattice boundaries and if so adjust
     * x or y lattice cell indices i */
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

    /* If the indices are outside the bound of the lattice */
    if (lat_x < 0 || lat_x >= _num_x ||
	lat_y < 0 || lat_y >= _num_y) {
        return NULL;
    }

    /* Compute local position of particle in the next level universe */
    double nextX = coords->getX() - (_origin.getX()
				     + (lat_x + 0.5) * _width_x);
    double nextY = coords->getY() - (_origin.getY()
				     + (lat_y + 0.5) * _width_y);

    /* Create a new localcoords object for the next level universe */
    LocalCoords* next_coords;

    if (coords->getNext() == NULL)
        next_coords = new LocalCoords(nextX, nextY);
    else
        next_coords = coords->getNext();

    short int universe_id = getUniverse(lat_x, lat_y)->getId();
    Universe* univ = universes.at(universe_id);
    next_coords->setUniverse(universe_id);

    /* Set lattice indices */
    coords->setLattice(_id);
    coords->setLatticeX(lat_x);
    coords->setLatticeY(lat_y);

    coords->setNext(next_coords);
    next_coords->setPrev(coords);

    /* Search the next lowest level universe for the cell */
    return univ->findCell(next_coords, universes);
}


/**
 * @brief Finds the next cell for a localcoords object along a trajectory 
 *        defined by some angle (in radians from 0 to PI). 
 * @details The method will update the localcoords passed in as an argument 
 *          to be the one at the boundary of the next cell crossed along the
 *          given trajectory. It will do this by recursively building a linked 
 *          list of localcoords from the localcoords passed in as an argument 
 *          down to the lowest level cell found. In the process it will set 
 *          the local coordinates for each localcoords in the linked list for 
 *          the lattice or universe that it is in. If the localcoords is 
 *          outside the bounds of the lattice or on the boundaries this method 
 *          will return NULL; otherwise it will return a pointer to the cell
 *          that the localcoords will reach next along its trajectory.
 * @param coords pointer to a localcoords object
 * @param angle the angle of the trajectory
 * @param universes a map of all of the universes passed in by the geometry
 * @return a pointer to a cell if found, NULL if no cell found
 */
Cell* Lattice::findNextLatticeCell(LocalCoords* coords, double angle,
				   std::map<short int, Universe*> universes) {

    /* Tests the upper, lower, left and right lattice cells adjacent to
     * the localcoord and uses the one with the shortest distance from
     * the current location of the localcoord */

    /* Initial distance is infinity */
    double distance = std::numeric_limits<double>::infinity();
    /* Current minimum distance */
    double d;

    /* Properties of the current localcoords */
    double x0 = coords->getX();
    double y0 = coords->getY();
    short int lattice_x = coords->getLatticeX();
    short int lattice_y = coords->getLatticeY();
    /* Slope of trajectory */
    double m = sin(angle) / cos(angle);

    /* Properties of the new location for localcoords */
    /* Current point of minimum distance */
    double x_curr, y_curr;
    /* x-coordinate on new lattice cell */
    double x_new = x0;
    /* y-coordinate on new lattice cell */
    double y_new = x0;
    /* New x lattice cell index */
    short int new_lattice_x;
    /* New y lattice cell index */
    short int new_lattice_y;
    /* Test point for computing distance */
    Point test;

    /* Check lower lattice cell Lower lattice cell */
    if (lattice_y >= 0 && angle >= M_PI) {
        y_curr = (lattice_y - _num_y/2.0) * _width_y;
	x_curr = x0 + (y_curr - y0) / m;
	test.setCoords(x_curr, y_curr);

	/* Check if the test point is within the bounds of the lattice */
	if (withinBounds(&test)) {
	    d = test.distanceToPoint(coords->getPoint());

	    /* Check if distance to test point is current minimum */
	    if (d < distance) {
	        distance = d;
		x_new = x_curr;
		y_new = y_curr;
	    }
	}
    }

    /* Upper lattice cell */
    if (lattice_y <= _num_y-1 && angle <= M_PI) {
        y_curr = (lattice_y - _num_y/2.0 + 1) * _width_y;
	x_curr = x0 + (y_curr - y0) / m;
	test.setCoords(x_curr, y_curr);

	/* Check if the test point is within the bounds of the lattice */
	if (withinBounds(&test)) {
	    d = test.distanceToPoint(coords->getPoint());

	    /* Check if distance to test point is current minimum */
	    if (d < distance) {
	        distance = d;
		x_new = x_curr;
		y_new = y_curr;
	    }
	}
    }

    /* Left lattice cell */
    if (lattice_x >= 0 && (angle >= M_PI/2 && angle <= 3*M_PI/2)) {
        x_curr = (lattice_x - _num_x/2.0) * _width_x;
	y_curr = y0 + m * (x_curr - x0);
	test.setCoords(x_curr, y_curr);

	/* Check if the test point is within the bounds of the lattice */
	if (withinBounds(&test)) {
	    d = test.distanceToPoint(coords->getPoint());

	    /* Check if distance to test point is current minimum */
	    if (d < distance) {
	        distance = d;
		x_new = x_curr;
		y_new = y_curr;
	    }
	}
    }

    /* Right lattice cell */
    if (lattice_x <= _num_x-1 && (angle <= M_PI/2 || angle >= 3*M_PI/2)) {
        x_curr = (lattice_x - _num_x/2.0 + 1) * _width_x;
	y_curr = y0 + m * (x_curr - x0);
	test.setCoords(x_curr, y_curr);
	
	/* Check if the test point is within the bounds of the lattice */
	if (withinBounds(&test)) {
	    d = test.distanceToPoint(coords->getPoint());

	    /* Check if distance to test point is current minimum */
	    if (d < distance) {
	        distance = d;
		x_new = x_curr;
		y_new = y_curr;
	    }
	}
    }

    /* If no point was found on the lattice cell, then the localcoords was
     * already on the boundary of the lattice */
    if (distance == INFINITY)
        return NULL;

    /* Otherwise a point was found inside a new lattice cell */
    else {
        /* Update the localcoords location to the point on the new lattice cell
	 * plus a small bit to ensure that its coordinates are inside cell */
        double delta_x = (x_new - coords->getX()) + cos(angle) * TINY_MOVE;
	double delta_y = (y_new - coords->getY()) + sin(angle) * TINY_MOVE;
	coords->adjustCoords(delta_x, delta_y);
	
	/* Compute the x and y indices for the new lattice cell */
	new_lattice_x = (int)floor((coords->getX() - _origin.getX())/_width_x);
	new_lattice_y = (int)floor((coords->getY() - _origin.getY())/_width_y);

	/* Check if the localcoord is on the lattice boundaries and if so adjust
	 * x or y lattice cell indices i */
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

	/* Check if new lattice cell indices are within the bounds, if not,
	 * new localcoords is now on the boundary of the lattice */
	if (new_lattice_x >= _num_x || new_lattice_x < 0)
	    return NULL;
	else if (new_lattice_y >= _num_y || new_lattice_y < 0)
	    return NULL;
	/* New localcoords is still within the interior of the lattice */
	else {
	    /* Update the localcoords lattice cell indices */
	    coords->setLatticeX(new_lattice_x);
	    coords->setLatticeY(new_lattice_y);

	    /* Move to next lowest level universe */
	    coords->prune();
	    Universe* univ = 
	        _universes.at(new_lattice_y).at(new_lattice_x).second;
	    LocalCoords* next_coords;

	    /* Compute local position of particle in next level universe */
	    double nextX = coords->getX() - (_origin.getX()
					   + (new_lattice_x + 0.5) * _width_x);
	    double nextY = coords->getY() - (_origin.getY()
					   + (new_lattice_y + 0.5) * _width_y);

	    /* Set the coordinates at the next level localcoord */
	    next_coords = new LocalCoords(nextX, nextY);
	    next_coords->setPrev(coords);
	    coords->setNext(next_coords);
	    
	    next_coords->setUniverse(univ->getId());
			
	    /* Search lower level universe */
	    return findCell(coords, universes);
	}
    }
}


/**
 * @brief Computes the flat source region base indices for each of the 
 *        lattice cells within this Lattice (ie, the minimum id for the flat 
 *        source regions within each lattice cell). Returns the number of 
 *        FSRs in the lattice.
 * @return the number of FSRs
 */
int Lattice::computeFSRMaps() {

    /* initialize a counter count */
    int count = 0;
    
    /* loop over universes in the lattice to set the map and update count */
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
 * @brief Converts a lattice's attributes to a character array representation.
 * @return character array of this lattice's attributes
 */
std::string Lattice::toString() {

    std::stringstream string;

    string << "Lattice id = " << _id << ", num cells along x = "
	   << _num_x << ", num cells along y = " << _num_y << ", x width = "
	   << _width_x << ", y width = " << _width_y;

    string << "\n\t\tUniverse ids within this lattice:\n\t\t";
    for (int i = _num_y-1; i > -1;  i--) {
        for (int j = 0; j < _num_x; j++)
	    string << _universes.at(i).at(j).first << "  ";
	string << "\n\t\t";
    }

    return string.str().c_str();
}


/**
 * @brief Prints a string representation of all of the lattice's objects to
 *        the console.
 */
void Lattice::printString() {
    log_printf(RESULT, toString().c_str());
}
