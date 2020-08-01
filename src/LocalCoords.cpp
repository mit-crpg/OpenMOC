#include "LocalCoords.h"

/**
 * @brief Constructor sets the x, y and z coordinates and position as a coord.
 * @param x the x-coordinate
 * @param y the y-coordinate
 * @param z the z-coordinate
 * @param first whether the LocalCoords is the first one, that will contain an
 *        array to all the next LocalCoords
 */
LocalCoords::LocalCoords(double x, double y, double z, bool first) {
  _coords.setCoords(x, y, z);
  _phi = 0.;
  _polar = M_PI_2;
  _universe = NULL;
  _lattice = NULL;
  _cell = NULL;
  _next = NULL;
  _prev = NULL;
  _version_num = 0;
  if (first) {
    _array_size = LOCAL_COORDS_LEN;
    _next_array = new LocalCoords[LOCAL_COORDS_LEN];
  }
  else {
    _array_size = 0;
    _next_array = NULL;
  }
  _position = -1;
}


/**
 * @brief Destructor.
 */
LocalCoords::~LocalCoords() {
  prune();
  if (_position == -1)
    deleteArray();
}


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
 * @brief Returns the direction angle in radians with respect to the z-axis.
 * @return the direction angle in radians
 */
double LocalCoords::getPolar() const {
  return _polar;
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
 * @brief Creates and returns a pointer to the next LocalCoords (nested
 *        deeper).
 * @param x the x-coordinate
 * @param y the y-coordinate
 * @param z the z-coordinate
 * @return pointer to the next LocalCoords that was just created
 */
LocalCoords* LocalCoords::getNextCreate(double x, double y, double z) {

  if (_next == NULL) {

    /* If more LocalCoords than the array size defined in constants.h,
    create a new one */
    if (_position + 1 >= _array_size) {
      _next = new LocalCoords(x, y, z, true);
      _next->setPrev(this);
    }
    else {
      _next = &_next_array[_position+1];
      _next->setPrev(this);
      _next->setNext(NULL);
      _next->setArrayPosition(_next_array, _position+1, _array_size);
      _next->getPoint()->setCoords(x, y, z);
      _next->setLattice(NULL);
      _next->setUniverse(NULL);
      _next->setCell(NULL);
      _next->setVersionNum(0);
    }
  }

  return _next;
}


/**
 * @brief Returns the LocalCoords position in the _next_array.
 * @return The position of this object in the underlying _next_array
 */
int LocalCoords::getPosition() {
  return _position;
}


/**
 * @brief Searches through the LocalCoords object to detect a loop.
 * @details A loop is assumed if the LocalCoords apparent length is greater
 *          1000 members
 */
void LocalCoords::detectLoop() {
  int n = 0;

  LocalCoords* iter = _next;
  while (iter != NULL) {
      iter = iter->getNext();
      n++;
      if (n > 1000)
        log_printf(ERROR, "Infinite loop of coords");
  }
  log_printf(DEBUG, "The LocalCoords is: %s\n", toString().c_str());
  log_printf(DEBUG, "The depth of the chain is %d \n", n);
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
 * @brief Returns the version of the LocalCoords object.
 * @details The version number differentiates otherwise matching FSR keys
 * @return The version number
 */
int LocalCoords::getVersionNum() {
  return _version_num;
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
 * @brief Sets the row index for the Lattice cell within which this
 *        LocalCoords resides.
 * @param lattice_x the row Lattice cell index
 */
void LocalCoords::setLatticeX(int lattice_x) {
  _lattice_x = lattice_x;
}


/**
 * @brief Sets the column index for the Lattice cell within which this
 *        LocalCoords resides.
 * @param lattice_y the column Lattice cell index
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
 * @brief Set the z-coordinate for this LocalCoords.
 * @param z the z-coordinate
 */
void LocalCoords::setZ(double z) {
  _coords.setZ(z);
}


/**
 * @brief Set the azimuthal angle for this LocalCoords.
 * @param phi the azimuthal angle, should be between 0 and 2 Pi
 */
void LocalCoords::setPhi(double phi) {
  _phi = phi;
}


/**
 * @brief Set the polar angle for this LocalCoords.
 * @param polar the polar angle
 */
void LocalCoords::setPolar(double polar) {
  _polar = polar;
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
 * @brief Sets the position of this LocalCoords in the array of LocalCoords,
 *        the pointer to this array (next LocalCoords for this) is also
 *        transfered
 * @param array pointer to the array of next LocalCoords
 * @param position index of this LocalCoords in said array
 * @param array_size size of the array
 */
void LocalCoords::setArrayPosition(LocalCoords* array, int position,
                                   int array_size) {
  _next_array = array;
  _position = position;
  _array_size = array_size;
}


/**
 * @brief Sets the version of the LocalCoords object.
 * @details The version number differentiates otherwise matching FSR keys
 * @param version_num the version number
 */
void LocalCoords::setVersionNum(int version_num) {
  _version_num = version_num;
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
 * @brief Translate all of the x,y,z coordinates for each LocalCoords object in
 *        the linked list.
 * @details This method will traverse the entire linked list and apply the
 *          translation to each element.
 * @param delta_x amount we wish to move x by
 * @param delta_y amount we wish to move y by
 * @param delta_z amount we wish to move z by
 */
void LocalCoords::adjustCoords(double delta_x, double delta_y, double delta_z) {

  /* Forward direction along linked list */
  LocalCoords* curr = this;
  while (curr != NULL) {
    curr->setX(curr->getX() + delta_x);
    curr->setY(curr->getY() + delta_y);
    curr->setZ(curr->getZ() + delta_z);
    curr = curr->getNext();
  }

  /* Reverse direction along linked list */
  curr = _prev;
  while (curr != NULL) {
    curr->setX(curr->getX() + delta_x);
    curr->setY(curr->getY() + delta_y);
    curr->setZ(curr->getZ() + delta_z);
    curr = curr->getPrev();
  }
  return;
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
  double delta_z = point->getZ() - curr->getZ();
  adjustCoords(delta_x, delta_y, delta_z);

  return;
}


/**
 * @brief Removes and frees memory for all LocalCoords beyond this one
 *        in the linked list.
 */
void LocalCoords::prune() {

  LocalCoords* curr = getLowestLevel();
  LocalCoords* next;

  /* Iterate over LocalCoords beneath this one in the linked list */
  while (curr != this) {
    next = curr->getPrev();
    if (curr->getPosition() == -1)
      curr->deleteArray();
    curr = next;
  }

  /* Set the next LocalCoord in the linked list to null */
  setNext(NULL);
}


/**
 * @brief Deletes the underlying array for next coordinates.
 */
void LocalCoords::deleteArray() {
  if (_next_array != NULL) {
      delete [] _next_array;
      _next_array = NULL;
  }
}


/**
 * @brief Copies a LocalCoords' values to this one.
 * @details Given a pointer to a LocalCoords, it first prunes it and then creates
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

    if (curr1 != NULL) {
      curr2 = curr2->getNextCreate(0, 0, 0);
    }
  }
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

    int univ_id = -1;
    if (curr->getUniverse() != NULL)
      univ_id = curr->getUniverse()->getId();
    int lat_id = -1;
    if (curr->getLattice() != NULL)
      lat_id = curr->getLattice()->getId();
    int cell_id = -1;
    if (curr->getCell() != NULL)
      cell_id = curr->getCell()->getId();

    if (curr->getType() == UNIV) {
      string << " UNIVERSE, x = " << curr->getX()
             << ", y = " << curr->getY()
             << ", z = " << curr->getZ();
      string << ", universe = " << univ_id;
      string << ", cell = " << cell_id;
    }
    else if (curr->getType() == LAT) {
      string << " LATTICE, x = " << curr->getX()
             << ", y = " << curr->getY()
             << ", z = " << curr->getZ()
             << ", universe = " << univ_id
             << ", lattice = " << lat_id
             << ", lattice_x = " << curr->getLatticeX()
             << ", lattice_y = " << curr->getLatticeY()
             << ", lattice_z = " << curr->getLatticeZ();
    }
    else {
      string << " NONE, x = " << curr->getX()
             << ", y = " << curr->getY()
             << ", z = " << curr->getZ()
             << ", universe = " << univ_id
             << ", lattice = " << lat_id
             << ", lattice_x = " << curr->getLatticeX()
             << ", lattice_y = " << curr->getLatticeY()
             << ", lattice_z = " << curr->getLatticeZ();
        string << ", cell = " << cell_id;
    }

    string << ", next:\n";
    curr = curr->getNext();
  }

  return string.str();
}
