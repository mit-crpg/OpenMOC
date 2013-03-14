/*
 * LocalCoords.cpp
 *
 *  Created on: Jan 25, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */


#include "LocalCoords.h"

/**
 * LocalCoords constructor
 * @param x the x-coordinate
 * @param y the y-coordinate
 */
LocalCoords::LocalCoords(double x, double y) {
	_coords.setCoords(x, y);
	_next = NULL;
	_prev = NULL;
}


/**
 * LocalCoords destructor recursively deletes all LocalCoords objects
 * behind this one in the linked list
 */
LocalCoords::~LocalCoords() { }


/**
 * Return the level (UNIV or LAT) of this localcoords
 * @return the level (UNIV or LAT)
 */
coordType LocalCoords::getType() {
	return _type;
}


/**
 * Return the universe id that this coordinate is in
 * @return the universe id
 */
short int LocalCoords::getUniverse() const {
    return _universe;
}


/**
 * Return the cell id for this localcoords
 * @return the cell id
 */
short int LocalCoords::getCell() const {
	return _cell;
}


/**
 * Return the lattice id that this coordinate is in
 * @return the lattice id
 */
short int LocalCoords::getLattice() const {
    return _lattice;
}

/**
 * Return lattice cell along the x-axis that this coordinate is in
 * @return lattice cell x
 */
short int LocalCoords::getLatticeX() const {
    return _lattice_x;
}


/**
 * Return lattice cell along the y-axis that this coordinate is in
 * @return lattice cell y
 */
short int LocalCoords::getLatticeY() const {
    return _lattice_y;
}


/**
 * Returns the x-coordinate
 * @return the x-coordinate
 */
double LocalCoords::getX() const {
    return _coords.getX();
}


/**
 * Returns the y-coordinate
 * @return the y-coordinate
 */
double LocalCoords::getY() const {
    return _coords.getY();
}


/**
 * Returns a pointer to the coordinates for this localcoord
 * @return pointer the coordinates
 */
Point* LocalCoords::getPoint() {
	return &_coords;
}


/**
 * Return a pointer to the localcoord at the next level if one exists
 * @return pointer to the next localcoord
 */
LocalCoords* LocalCoords::getNext() const {
    return _next;
}


LocalCoords* LocalCoords::getPrev() const {
	return _prev;
}


/**
 * Set the level for this localcoords
 * @param level the level for this localcoords (UNIV or LAT)
 */
void LocalCoords::setType(coordType type) {
	_type = type;
}


/**
 * Sets the universe id that this coordinate is in
 * @param universe the universe id
 */
void LocalCoords::setUniverse(short int universe) {
    _universe = universe;
}


/**
 * Set the cell id for this localcoords
 * @param cell the cell id
 */
void LocalCoords::setCell(short int cell) {
	_cell = cell;
}


/**
 * Sets the lattice id that this coordinate is in
 */
void LocalCoords::setLattice(short int lattice) {
    _lattice = lattice;
}


/**
 * Sets the lattice cell along the x-axis that this coordinate is in
 * @param lattice_x the x lattice cell
 */
void LocalCoords::setLatticeX(short int lattice_x) {
    _lattice_x = lattice_x;
}


/**
 * Sets the lattice cell along the y-axis that this coordinate is in
 * @param lattice_y the y lattice cell
 */
void LocalCoords::setLatticeY(short int lattice_y) {
    _lattice_y = lattice_y;
}


/**
 * Set the x-coordinate
 * @param x the x-coordinate
 */
void LocalCoords::setX(double x) {
	_coords.setX(x);
}


/**
 * Set the y-coordinate
 * @param y the y-coordinate
 */
void LocalCoords::setY(double y) {
	_coords.setY(y);
}


/**
 * Sets the pointer to the next localcoord in the linked list
 * @param next pointer to the next localcoord
 */
void LocalCoords::setNext(LocalCoords* next) {
    _next = next;
}


/**
 * Sets the pointer to the previous localcoord in the linked list
 * @param prev pointer to the previous localcoord
 */
void LocalCoords::setPrev(LocalCoords* prev) {
	_prev = prev;
}


/**
 * Find and return the last localcoord in the linked list wich represents
 * the local coordinates on the lowest level of a geometry of nested universes
 * @return a pointer to the last localcoord object in the list
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
 * Update all of the x,y coordinates for each localcoord object in the linked
 * list. This method will traverse the entire linked list and apply the
 * translation to each element
 * @param delta_x amount we wish to move x by
 * @param delta y amount we wish to move y by
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
 * Update the last element in the linked list (the one at the lowest level
 * of nested universes) to have the same coordinates as a given point
 * @param point a pointer to a point
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
 * Removes and frees memory for all localcoords beyond this one
 * in the linked list
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
 * Given a pointer to a localcoords, it first prunes it and then creates
 * a copy of the linked list of localcoords in the linked list below this one
 * to give to the input localcoords
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
 * Converts this localcoords's attributes to a character array representation
 * @return a character array of its member's attributes
 */
std::string LocalCoords::toString() {

	std::stringstream string;
	LocalCoords* curr = this;

	/* Loops over all localcoords lower than this one in the list */
	while (curr != NULL) {
		string << "LocalCoords: level = ";

		if (curr->getType() == UNIV) {
			string << " UNIVERSE, x = " << curr->getX() << ", y = " << curr->getY()
					<< ", universe = " << curr->getUniverse() << ", cell = " <<
					curr->getCell();
		}
		else if (curr->getType() == LAT){
			string << " LATTICE, x = " << curr->getX() << ", y = " << curr->getY()
				<< ", universe = " << curr->getUniverse() << ", lattice = " <<
				curr->getLattice() << ", lattice_x = " << curr->getLatticeX()
				<< ", lattice_y = " << curr->getLatticeY();
		}
		else {
			string << " NONE, x = " << curr->getX() << ", y = " << curr->getY()
				<< ", universe = " << curr->getUniverse() << ", lattice = " <<
				curr->getLattice() << ", lattice_x = " << curr->getLatticeX()
				<< ", lattice_y = " << curr->getLatticeY()
				<< ", cell = " << curr->getCell();
		}

		string << ", next:\n";
		curr = curr->getNext();
	}

	return string.str();
}
