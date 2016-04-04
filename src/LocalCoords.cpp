#include "LocalCoords.h"

/**
 * @brief Constructor sets the x and y coordinates.
 * @param x the x-coordinate
 * @param y the y-coordinate
 */
LocalCoords::LocalCoords(double x, double y, double z) {
  _coords.setCoords(x, y, z);
  _phi = 0.;
  _universe = NULL;
  _lattice = NULL;
  _cell = NULL;
  _next = NULL;
  _prev = NULL;
}


/**
 * @brief Destructor.
 */
LocalCoords::~LocalCoords() { }


/**
 * @brief Return the level (UNIV or LAT) of this LocalCoords.
 * @return the nested Universe level (UNIV or LAT)
 */
coordType LocalCoords::getType() {
  return _type;
}


/**
 * @brief Return the Universe within which this LocalCoords resides.
 * @return the Universe
 */
Universe* LocalCoords::getUniverse() const {
  return _universe;
}


/**
 * @brief Return the Cell within which this LocalCoords resides.
 * @return the Cell
 */
Cell* LocalCoords::getCell() const {
  return _cell;
}


/**
 * @brief Return the Lattice within which this LocalCoords resides.
 * @return the Lattice
 */
Lattice* LocalCoords::getLattice() const {
  return _lattice;
}

/**
 * @brief Return the first index of the Lattice cell within which this
 *        LocalCoords resides.
 * @return the first Lattice cell index
 */
int LocalCoords::getLatticeX() const {
  return _lattice_x;
}


/**
 * @brief Return the second index of the Lattice cell within which this
 *        LocalCoords resides.
 * @return the second Lattice cell index
 */
int LocalCoords::getLatticeY() const {
  return _lattice_y;
}


/**
 * @brief Return the third index of the Lattice cell within which this
 *        LocalCoords resides.
 * @return the third Lattice cell index
 */
int LocalCoords::getLatticeZ() const {
  return _lattice_z;
}


/**
 * @brief Returns the x-coordinate for this LocalCoords location.
 * @return the x-coordinate of this LocalCoords location
 */
double LocalCoords::getX() const {
  return _coords.getX();
}


/**
 * @brief Returns the y-coordinate for this LocalCoords location.
 * @return the y-coordinate of this LocalCoords location
 */
double LocalCoords::getY() const {
  return _coords.getY();
}


/**
 * @brief Returns the z-coordinate for this LocalCoords location.
 * @return the z-coordinate of this LocalCoords location
 */
double LocalCoords::getZ() const {
  return _coords.getZ();
}


/**
 * @brief Returns the direction angle in radians with respect to the x-axis.
 * @return the direction angle in radians
 */
double LocalCoords::getPhi() const {
  return _phi;
}


/**
 * @brief Returns a pointer to the Point containing the coordinates for this
 *        LocalCoord.
 * @return pointer to the Point containing the x and y coordinates
 */
Point* LocalCoords::getPoint() {
  return &_coords;
}


/**
 * @brief Return a pointer to the LocalCoord at the next lower nested Universe
 *        level if one exists.
 * @return pointer to the next LocalCoord
 */
LocalCoords* LocalCoords::getNext() const {
  return _next;
}


/**
 * @brief Return a pointer to the LocalCoord at the next higher nested Universe
 *        level if one exists.
 * @return pointer to the previous LocalCoord
 */
LocalCoords* LocalCoords::getPrev() const {
  return _prev;
}


/**
 * @brief Set the type of LocalCoords (UNIV or LAT).
 * @param type the type for LocalCoords (UNIV or LAT)
 */
void LocalCoords::setType(coordType type) {
  _type = type;
}


/**
 * @brief Set the Universe within which this LocalCoords resides.
 * @param universe the Universe
 */
void LocalCoords::setUniverse(Universe* universe) {
  _universe = universe;
}


/**
 * @brief Set the Cell within which this LocalCoords resides.
 * @param cell the Cell
 */
void LocalCoords::setCell(Cell* cell) {
  _cell = cell;
}


/**
 * @brief Sets the Lattice within which this LocalCoords resides.
 * @param lattice the Lattice
 */
void LocalCoords::setLattice(Lattice* lattice) {
  _lattice = lattice;
}


/**
 * @brief Sets the x index for the Lattice cell within which this
 *        LocalCoords resides.
 * @param lattice_x the x Lattice cell index
 */
void LocalCoords::setLatticeX(int lattice_x) {
  _lattice_x = lattice_x;
}


/**
 * @brief Sets the y index for the Lattice cell within which this
 *        LocalCoords resides.
 * @param lattice_y the y Lattice cell index
 */
void LocalCoords::setLatticeY(int lattice_y) {
  _lattice_y = lattice_y;
}


/**
 * @brief Sets the z index for the Lattice cell within which this
 *        LocalCoords resides.
 * @param lattice_z the z Lattice cell index
 */
void LocalCoords::setLatticeZ(int lattice_z) {
  _lattice_z = lattice_z;
}


/**
 * @brief Set the x-coordinate for this LocalCoords.
 * @param x the x-coordinate
 */
void LocalCoords::setX(double x) {
  _coords.setX(x);
}


/**
 * @brief Set the y-coordinate for this Localcoords.
 * @param y the y-coordinate
 */
void LocalCoords::setY(double y) {
  _coords.setY(y);
}


/**
 * @brief Set the z-coordinate for this Localcoords.
 * @param z the z-coordinate
 */
void LocalCoords::setZ(double z) {
  _coords.setZ(z);
}


/**
 * @brief Set the direction angle in radians for this LocalCoords.
 * @param angle the direction angle in radians
 */
void LocalCoords::setPhi(double phi) {
  _phi = phi;
}


/**
 * @brief Increment the direction angle in radians for this LocalCoords.
 * @param phi the incremental direction angle in radians
 */
void LocalCoords::incrementPhi(double phi) {
  _phi += phi;
}


/**
 * @brief Sets the pointer to the LocalCoords on the next lower nested Universe
 *        level.
 * @param next pointer to the next LocalCoords
 */
void LocalCoords::setNext(LocalCoords* next) {
  _next = next;
}


/**
 * @brief Sets the pointer to the LocalCoords on the next higher nested
 *        Universe level.
 * @param prev pointer to the previous LocalCoords
 */
void LocalCoords::setPrev(LocalCoords* prev) {
  _prev = prev;
}


/**
 * @brief Find and return the last LocalCoords in the linked list which
 *        represents the local coordinates on the lowest level of a geometry
 *        of nested universes.
 * @details Traverses a linked list of LocalCoords to find the one at the
 *          lowest nested Universe level.
 * @return a pointer to the last LocalCoords object in the list
 */
LocalCoords* LocalCoords::getLowestLevel() {
  LocalCoords* curr = this;

  /* Traverse linked list */
  while (curr->getNext() != NULL)
    curr = curr->getNext();

  return curr;
}


/**
 * @brief Find and return the first LocalCoords in the linked list which
 *        represents the local coordinates on the highest level of a geometry
 *        of nested universes.
 * @details Traverses a linked list of LocalCoords to find the one at the
 *          highest nested Universe level.
 * @return a pointer to the first LocalCoords object in the list
 */
LocalCoords* LocalCoords::getHighestLevel() {
  LocalCoords* curr = this;

  /* Traverse linked list */
  while (curr->getPrev() != NULL)
    curr = curr->getPrev();

  return curr;
}


/**
 * @brief Translate all of the x,y coordinates for each LocalCoords object in
 *        the linked list.
 * @details This method will traverse the entire linked list and apply the
 *          translation to each element.
 * @param delta amount we wish to move the point by
 */
void LocalCoords::adjustCoords(double delta) {
  double new_x = getX() + cos(_phi) * delta;
  double new_y = getY() + sin(_phi) * delta;
  setX(new_x);
  setY(new_y);

  if (_next != NULL)
    _next->adjustCoords(delta);
}


/**
 * @brief Update the last element in the linked list (the one at the lowest
 *        level of nested Universes) to have the same coordinates as a
 *        given point.
 * @param point a pointer to a point of interest
 */
void LocalCoords::updateMostLocal(Point* point) {

  /* Get the lowest level coordinate */
  LocalCoords* curr = getLowestLevel();

  /* Translate coordinates by appropriate amount */
  double delta_x = point->getX() - curr->getX();
  double delta_y = point->getY() - curr->getY();
  double delta = sqrt(delta_x*delta_x + delta_y*delta_y);
  adjustCoords(delta);

  return;
}


/**
 * @brief Removes and frees memory for all LocalCoords beyond this one
 *        in the linked list
 */
void LocalCoords::prune() {

  LocalCoords* curr = getLowestLevel();
  LocalCoords* next = curr->getPrev();

  /* Iterate over LocalCoords beneath this one in the linked list */
  while (curr != this) {
    next = curr->getPrev();
    delete curr;
    curr = next;
  }

  /* Set the next LocalCoord in the linked list to null */
  setNext(NULL);
}


/**
 * @brief Copies a LocalCoords' values to this one.
 * details Given a pointer to a LocalCoords, it first prunes it and then creates
 *         a copy of the linked list of LocalCoords in the linked list below
 *         this one to give to the input LocalCoords.
 * @param coords a pointer to the LocalCoords to give the linked list copy to
 */
void LocalCoords::copyCoords(LocalCoords* coords) {

  LocalCoords* curr1 = this;
  LocalCoords* curr2 = coords;

  /* Prune the LocalCoords linked list */
  curr2->prune();

  /* Iterate over this LocalCoords linked list and create a
   * copy of it for the input LocalCoords */
  while (curr1 != NULL) {
    curr2->setX(curr1->getX());
    curr2->setY(curr1->getY());
    curr2->setZ(curr1->getZ());
    curr2->setPhi(curr1->getPhi());
    curr2->setUniverse(curr1->getUniverse());

    if (curr1->getType() == UNIV) {
      curr2->setType(UNIV);
      curr2->setCell(curr1->getCell());
    }
    else {
      curr2->setLattice(curr1->getLattice());
      curr2->setLatticeX(curr1->getLatticeX());
      curr2->setLatticeY(curr1->getLatticeY());
      curr2->setLatticeZ(curr1->getLatticeZ());
      curr2->setType(LAT);
    }

    curr1 = curr1->getNext();

    if (curr1 != NULL && curr2->getNext() == NULL) {
      LocalCoords* new_coords = new LocalCoords(0.0, 0.0, 0.0);
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
 * @brief Converts this LocalCoords's attributes to a character array
 *        representation.
 * @return a character array of the LocalCoord's attributes
 */
std::string LocalCoords::toString() {

  std::stringstream string;
  LocalCoords* curr = this;

  /* Loops over all LocalCoords lower than this one in the list */
  while (curr != NULL) {
    string << "LocalCoords: level = ";

    if (curr->getType() == UNIV) {
      string << " UNIVERSE, x = " << curr->getX()
             << ", y = " << curr->getY()
             << ", z = " << curr->getZ()
             << ", universe = " << curr->getUniverse()->getId()
             << ", cell = " << curr->getCell()->getId();
    }
    else if (curr->getType() == LAT) {
      string << " LATTICE, x = " << curr->getX()
             << ", y = " << curr->getY()
             << ", z = " << curr->getZ()
             << ", universe = " << curr->getUniverse()->getId()
             << ", lattice = " << curr->getLattice()->getId()
             << ", lattice_x = " << curr->getLatticeX()
             << ", lattice_y = " << curr->getLatticeY()
             << ", lattice_z = " << curr->getLatticeZ();
    }
    else {
      string << " NONE, x = " << curr->getX()
             << ", y = " << curr->getY()
             << ", z = " << curr->getZ()
             << ", universe = " << curr->getUniverse()->getId()
             << ", lattice = " << curr->getLattice()->getId()
             << ", lattice_x = " << curr->getLatticeX()
             << ", lattice_y = " << curr->getLatticeY()
             << ", lattice_z = " << curr->getLatticeZ()
             << ", cell = " << curr->getCell();
    }

    string << ", next:\n";
    curr = curr->getNext();
  }

  return string.str();
}
