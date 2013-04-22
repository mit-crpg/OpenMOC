#include "LocalCoords.h"

/**
 * @brief Constructor sets the x and y coordinates.
 * @param x the x-coordinate
 * @param y the y-coordinate
 */
LocalCoords::LocalCoords(double x, double y) {
    _coords.setCoords(x, y);
    _next = NULL;
    _prev = NULL;
}


/**
 * @brief Destructor.
 */
LocalCoords::~LocalCoords() { }


/**
 * @brief Return the level (UNIV or LAT) of this localcoords.
 * @return the level (UNIV or LAT)
 */
coordType LocalCoords::getType() {
    return _type;
}


/**
 * @brief Return the id of the universe within which this localcoords resides.
 * @return the universe id
 */
short int LocalCoords::getUniverse() const {
    return _universe;
}


/**
 * @brief Return the id of the cell within which this localcoords resides.
 * @return the cell id
 */
short int LocalCoords::getCell() const {
    return _cell;
}


/**
 * @brief Return the id of the lattice within which this localcoords resides.
 * @return the lattice id
 */
short int LocalCoords::getLattice() const {
    return _lattice;
}

/**
 * @brief Return the first index of the lattice cell within which this 
 *        localcoords resides.
 * @return the first lattice cell index
 */
short int LocalCoords::getLatticeX() const {
    return _lattice_x;
}


/**
 * @brief Return the second index of the lattice cell within which this
 *        localcoords resides.
 * @return the second lattice cell index
 */
short int LocalCoords::getLatticeY() const {
    return _lattice_y;
}


/**
 * @brief Returns the x-coordinate for this localcoords location.
 * @return the x-coordinate of this localcoords location
 */
double LocalCoords::getX() const {
    return _coords.getX();
}


/**
 * @brief Returns the y-coordinate for this localcoords location.
 * @return the y-coordinate of this localcoords location
 */
double LocalCoords::getY() const {
    return _coords.getY();
}


/**
 * @brief Returns a pointer to the point containing the coordinates for this
 *        localcoord.
 * @return pointer to the point containing the x and y coordinates
 */
Point* LocalCoords::getPoint() {
    return &_coords;
}


/**
 * @brief Return a pointer to the localcoord at the next lower nested universe
 *        level if one exists.
 * @return pointer to the next localcoord
 */
LocalCoords* LocalCoords::getNext() const {
    return _next;
}


/**
 * @brief Return a pointer to the localcoord at the next higher nested universe
 *        level if one exists.
 * @return pointer to the previous localcoord
 */
LocalCoords* LocalCoords::getPrev() const {
    return _prev;
}


/**
 * @brief Set the type of localcoords (UNIV or LAT).
 * @param type the type for localcoords (UNIV or LAT)
 */
void LocalCoords::setType(coordType type) {
    _type = type;
}


/**
 * @brief Set the id of the universe within which this localcoords resides.
 * @param universe the universe id
 */
void LocalCoords::setUniverse(short int universe) {
    _universe = universe;
}


/**
 * @brief Set the id of the cell within which this localcoords resides.
 * @param cell the cell id
 */
void LocalCoords::setCell(short int cell) {
    _cell = cell;
}


/**
 * @brief Sets the id of the lattice within which this localcoords resides.
 * @param lattice the lattice id
 */
void LocalCoords::setLattice(short int lattice) {
    _lattice = lattice;
}


/**
 * @brief Sets the first index for the lattice cell within which this 
 *        localcoords resides.
 * @param lattice_x the first lattice cell index
 */
void LocalCoords::setLatticeX(short int lattice_x) {
    _lattice_x = lattice_x;
}


/**
 * @brief Sets the second index for the lattice cell within which this 
 *        localcoords resides.
 * @param lattice_y the second lattice cell index
 */
void LocalCoords::setLatticeY(short int lattice_y) {
    _lattice_y = lattice_y;
}


/**
 * @brief Set the x-coordinate for this localcoords.
 * @param x the x-coordinate
 */
void LocalCoords::setX(double x) {
    _coords.setX(x);
}


/**
 * @brief Set the y-coordinate for this localcoords.
 * @param y the y-coordinate
 */
void LocalCoords::setY(double y) {
    _coords.setY(y);
}


/**
 * @brief Sets the pointer to the localcoords on the next lower nested universe
 *        level.
 * @param next pointer to the next localcoords
 */
void LocalCoords::setNext(LocalCoords* next) {
    _next = next;
}


/**
 * @brief Sets the pointer to the localcoords on the next higher nested
 *        universe level.
 * @param prev pointer to the previous localcoords
 */
void LocalCoords::setPrev(LocalCoords* prev) {
    _prev = prev;
}


/**
 * @brief Find and return the last localcoords in the linked list wich 
 *        represents the local coordinates on the lowest level of a geometry
 *        of nested universes.
 * @details Traverses a linked list of localcoords to find the one at the
 *          lowest universe level.
 * @return a pointer to the last localcoords object in the list
 */
LocalCoords* LocalCoords::getLowestLevel() {
    LocalCoords* curr = this;

    if (curr)

    /* Traverse linked list */
    while (curr->getNext() != NULL)
        curr = curr->getNext();

    return curr;
}


/**
 * @brief Translate all of the x,y coordinates for each localcoords object in 
 *        the linked list. 
 * @details This method will traverse the entire linked list and apply the 
 *          translation to each element.
 * @param delta_x amount we wish to move x by
 * @param delta_y amount we wish to move y by
 */
void LocalCoords::adjustCoords(double delta_x, double delta_y) {

    /* Forward direction along linked list */
    LocalCoords* curr = this;
    while (curr != NULL) {
        curr->setX(curr->getX() + delta_x);
	curr->setY(curr->getY() + delta_y);
	curr = curr->getNext();
    }

    /* Reverse direction along linked list */
    curr = _prev;
    while (curr != NULL) {
        curr->setX(curr->getX() + delta_x);
	curr->setY(curr->getY() + delta_y);
	curr = curr->getPrev();
    }
    return;
}


/**
 * @brief Update the last element in the linked list (the one at the lowest 
 *        level of nested universes) to have the same coordinates as a 
 *        given point.
 * @param point a pointer to a point of interest
 */
void LocalCoords::updateMostLocal(Point* point) {

    /* Get the lowest level coordinate */
    LocalCoords* curr = getLowestLevel();

    /* Translate coordinates by appropriate amount */
    double delta_x = point->getX() - curr->getX();
    double delta_y = point->getY() - curr->getY();
    adjustCoords(delta_x, delta_y);

    return;
}


/**
 * @brief Removes and frees memory for all localcoords beyond this one
 *        in the linked list
 */
void LocalCoords::prune() {

    LocalCoords* curr = getLowestLevel();
    LocalCoords* next = curr->getPrev();

    /* Iterate over localcoords beneath this one in the linked list */
    while (curr != this) {
        next = curr->getPrev();
	delete curr;
	curr = next;
    }

    /* Set the next localcoord in the linked list to null */
    setNext(NULL);
}


/**
 * @brief Copies a localcoords' values to this one.
 * details Given a pointer to a localcoords, it first prunes it and then creates
 *         a copy of the linked list of localcoords in the linked list below 
 *         this one to give to the input localcoords.
 * @param coords a pointer to the localcoords to give the linked list copy to
 */
void LocalCoords::copyCoords(LocalCoords* coords) {

    LocalCoords* curr1 = this;
    LocalCoords* curr2 = coords;

    /* Prune the localcoords linked list */
    curr2->prune();

    /* Iterate over this localcoords linked list and create a
     * copy of it for the input localcoords */
    while (curr1 != NULL) {
        curr2->setX(curr1->getX());
	curr2->setY(curr1->getY());
	curr2->setUniverse(curr1->getUniverse());

	if (curr1->getType() == UNIV) {
	    curr2->setType(UNIV);
	    curr2->setCell(curr1->getCell());
	}
	else {
	    curr2->setLattice(curr1->getLattice());
	    curr2->setLatticeX(curr1->getLatticeX());
	    curr2->setLatticeY(curr1->getLatticeY());
	    curr2->setType(LAT);
	}

	curr1 = curr1->getNext();
	
	if (curr1 != NULL && curr2->getNext() == NULL) {
	    LocalCoords* new_coords = new LocalCoords(0.0, 0.0);
	    curr2->setNext(new_coords);
	    new_coords->setPrev(curr2);
	    curr2 = new_coords;
	}
	else if (curr1 != NULL)
	    curr2 = curr2->getNext();
    }

    /* Prune any remainder from the old coords linked list */
    if (curr2 != NULL)
        curr2->prune();
}


/**
 * @brief Converts this localcoords's attributes to a character array 
 *        representation.
 * @return a character array of its member's attributes
 */
std::string LocalCoords::toString() {

    std::stringstream string;
    LocalCoords* curr = this;

    /* Loops over all localcoords lower than this one in the list */
    while (curr != NULL) {
        string << "LocalCoords: level = ";

	if (curr->getType() == UNIV) {
	    string << " UNIVERSE, x = " << curr->getX() << ", y = " 
		   << curr->getY() << ", universe = " << curr->getUniverse() 
		   << ", cell = " << curr->getCell();
	}
	else if (curr->getType() == LAT){
	    string << " LATTICE, x = " << curr->getX() << ", y = " 
		   << curr->getY() << ", universe = " << curr->getUniverse() 
		   << ", lattice = " << curr->getLattice() << ", lattice_x = " 
		   << curr->getLatticeX() << ", lattice_y = " 
		   << curr->getLatticeY();
	}
	else {
	    string << " NONE, x = " << curr->getX() << ", y = " << curr->getY()
		   << ", universe = " << curr->getUniverse() << ", lattice = " 
		   << curr->getLattice() << ", lattice_x = " 
		   << curr->getLatticeX() << ", lattice_y = " 
		   << curr->getLatticeY() << ", cell = " << curr->getCell();
	}

	string << ", next:\n";
        curr = curr->getNext();
    }

    return string.str();
}
