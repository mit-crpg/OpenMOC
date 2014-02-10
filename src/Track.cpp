#include "Track.h"


/*
 * @brief Constructor initializes an empty Track.
 */
Track::Track() { }



/**
 * @brief Destructor clears the Track segments container.
 */
Track::~Track() {
  clearSegments();
}


/**
 * @brief Set the values for the Track's start and end point and angle.
 * @param start_x the x-coordinate at the starting point
 * @param start_y the y-coordinate at the starting point
 * @param end_x the x-coordinate at the ending point
 * @param end_y the y-coordinate at the ending point
 * @param phi the track's azimuthal angle (\f$ \theta \in [0, \pi] \f$)
 */
 void Track::setValues(const double start_x, const double start_y,
                       const double end_x, const double end_y,
                       const double phi) {
   _start.setCoords(start_x, start_y);
   _end.setCoords(end_x, end_y);
   _phi = phi;
}


/**
 * @brief Initializes a Track's unique ID.
 * @details This is set by the trackgenerator to correspond to the Track's
 *          location in a 2D ragged array of all tracks.
 * @param uid the Track's unique ID
 */
void Track::setUid(int uid) {
  _uid = uid;
}

/**
 * @brief Set the Track's azimuthal angle.
 * @param phi the azimuthal angle
 */
void Track::setPhi(const double phi) {
  _phi = phi;
}


/**
 * @brief Set the index for the Track's azimuthal angle index.
 * @details The azimuthal angle index corresponds to a an array of all
 *          azimuthal angles for \f$ \theta \in [0, \pi] \f$ owned by
 *          the TrackGenerator class.
 * @param index the azimuthal angle index
 */
void Track::setAzimAngleIndex(const int index) {
  _azim_angle_index = index;
}


/**
 * @brief Adds a segment pointer to this Track's list of segments.
 * @details This method assumes that segments are added in order of their
 *          starting location from the Track's start point.
 * @param segment a pointer to the segment
 */
void Track::addSegment(segment* segment) {

  try {
    _segments.push_back(*segment);
  }
  catch (std::exception &e) {
      log_printf(ERROR, "Unable to add a segment to Track. Backtrace:"
                 "\n%s", e.what());
  }
}


/**
 * @brief Sets the direction in which the flux leaving this Track along its
 *        "forward" direction is passed to reflective Track for boundary
 *        conditions.
 * @details Sets whether or not to pass the outgoing flux from this Track
 *          along its "forward" direction to the "forward" direction (false)
 *          or "reverse" direction (true) of the Track reflecting out of this
 *          one at the boundary. This is used for reflective boundary
 *          conditions.
 * @param refl_in the "forward" (false) or "reverse (true) direction
 */
void Track::setReflIn(const bool refl_in) {
  _refl_in = refl_in;
}


/**
 * @brief Sets the direction in which the flux leaving this Track along
 *        its "reverse" direction is passed to reflective track for boundary
 *        conditions.
 * @details Sets whether or not to pass the outgoing flux from this Track
 *          along its "reverse" direction to the "forward" direction (false)
 *          or "reverse" direction (true) of the Track reflecting out of this
 *          one at the boundary. This is used for reflective boundary
 *          conditions.
 * @param refl_out "forward" (false) or "reverse (true) direction
 */
void Track::setReflOut(const bool refl_out) {
  _refl_out = refl_out;
}


/**
 * @brief Sets the boundary condition for the incoming flux along the Track's
 *        "forward" direction.
 * @details The boolean represents vacuum (false) or reflective (true)
 *          boundary conditions.
 * @param bc_in boundary condition for the incoming flux in the "forward"
 *        direction
 */
void Track::setBCIn(const bool bc_in) {
  _bc_in = bc_in;
}


/**
 * @brief Sets the boundary condition for the incoming flux along the Track's
 *        "reverse" direction.
 * @details The boolean represents vacuum (false) or reflective (true)
 *          boundary conditions.
 * @param bc_out boundary condition for the incoming flux in the "reverse"
 *        direction
 */
void Track::setBCOut(const bool bc_out) {
  _bc_out = bc_out;
}


/**
 * @brief Sets the track reflecting into this Track's "forward" direction.
 * @param track_in pointer to the Track reflecting into the "forward" direction
 */
void Track::setTrackIn(Track* track_in) {
  _track_in = track_in;
}


/**
 * @brief Sets the track reflecting into this Track's "reverse" direction.
 * @param track_out pointer to the Track reflecting into the "reverse" direction
 */
void Track::setTrackOut(Track* track_out) {
  _track_out = track_out;
}


/**
 * @brief Sets the first index of the track reflecting into this Track's
 *        "forward" direction in the 2D jagged array of tracks.
 * @param i the first index of the incoming track along the "forward" direction
 */
void Track::setTrackInI(int i) {
  _track_in_i = i;
}


/**
 * @brief Sets the second index of the track reflecting into this Track's
 *        "forward" direction in the 2D jagged array of tracks.
 * @param j the second index of the incoming track along the "forward" direction
 */
void Track::setTrackInJ(int j) {
  _track_in_j = j;
}


/**
 * @brief Sets the first index of the track reflecting into this Track's
 *        "reverse" direction in the 2D jagged array of Tracks.
 * @param i the first index of the incoming Track along the "reverse" direction
 */
void Track::setTrackOutI(int i) {
  _track_out_i = i;
}


/**
 * @brief Sets the second index of the Track reflecting into this Track's
 *        "reverse" direction in the 2D jagged array of tracks.
 * @param j the second index of the incoming Track along the "reverse" direction
 */
void Track::setTrackOutJ(int j) {
  _track_out_j = j;
}


/**
 * @brief Returns whether to give the outgoing flux to the "forward" (false) or
 *        "reverse" (true) direction of the Track reflecting out of this one
 *        along its "forward" direction.
 * @return "forward" (false) "reverse" (true) direction of outgoing Track
 */
bool Track::isReflIn() const {
  return _refl_in;
}


/**
 * @brief Returns whether to give the outgoing flux to the "forward" (false) or
 *        "reverse" (true) direction of the Track reflecting out of this one
 *        along its "reverse" direction.
 * @return "forward" (false) "reverse" (true) direction of outgoing Track
 */
bool Track::isReflOut() const {
  return _refl_out;
}


/**
 * @brief Returns the boundary condition for the flux along the Track's
 *        "forward" direction.
 * @return vacuum (false) or reflective (true) reflective boundary conditions
 */
bool Track::getBCIn() const {
  return _bc_in;
}


/**
 * @brief Returns the boundary condition for the flux along the Track's
 *        "reverse" direction.
 * @return vacuum (false) or reflective (true) reflective boundary conditions
 */
bool Track::getBCOut() const {
  return _bc_out;
}


/**
 * @brief Returns a pointer to the Track's end Point.
 * @return a pointer to the Track's end Point
 */
Point* Track::getEnd() {
  return &_end;
}


/**
 * @brief Returns a pointer to the Track's start Point.
 * @return a pointer to the Track's start Point
 */
Point* Track::getStart() {
  return &_start;
}


/**
 * @brief Return the Track's azimuthal angle (with respect to the x-axis).
 * @return the azimuthal angle \f$ \theta \in [0, \pi] \f$
 */
double Track::getPhi() const {
  return _phi;
}


/**
 * @brief Return the index for the Track's azimuthal angle (with respect to the
 *        x-axis).
 * @return th azimuthal angle index
 */
int Track::getAzimAngleIndex() const {
  return _azim_angle_index;
}


/**
 * @brief Returns the first index of the Track reflecting out of this one along
 *        its "forward" direction in the 2D jagged array of all Tracks.
 * @return the first index of the reflecting Track
 */
int Track::getTrackInI() const {
  return _track_in_i;
}


/**
 * @brief Returns the second index of the Track reflecting out of this one along
 *        its "forward" direction in the 2D jagged array of all Tracks.
 * @return the second index of the reflecting Track
 */
int Track::getTrackInJ() const {
  return _track_in_j;
}



/**
 * @brief Returns the first index of the Track reflecting out of this one along
 *        its "reverse" direction in the 2D jagged array of all Tracks.
 * @return the first index of the reflecting Track
 */
int Track::getTrackOutI() const {
  return _track_out_i;
}


/**
 * @brief Returns the second index of the Track reflecting out of this one along
 *        its "reverse" direction in the 2D jagged array of all Tracks.
 * @return the second index of the reflecting Track
 */
int Track::getTrackOutJ() const {
  return _track_out_j;
}


/**
 * @brief Checks whether a Point is contained along this Track.
 * @param point a pointer to the Point of interest
 * @return true if the Point is on the Track, false otherwise
 */
bool Track::contains(Point* point) {

  /* The slope of the Track */
  double m;

  /* The distance from the Point to the Track */
  double dist;

  /* If the Point is outside of the bounds of the start and end Points of the
   * Track it does not lie on the Track */
  if (!(((point->getX() <= _start.getX()+1.0E-2 &&
          point->getX() >= _end.getX()-1.0E-2)
      || (point->getX() >= _start.getX()-1.0E-2 &&
          point->getX() <= _end.getX()+1.0E-2)) &&
        ((point->getY() <= _start.getY()+1.0E-2 &&
          point->getY() >= _end.getY()-1.0E-2)
      || (point->getY() >= _start.getY()-1.0E-2 &&
          point->getY() <= _end.getY()+1.0E-2)))) {

    return false;
  }


  /* If the Track is vertical */
  if (fabs(_phi - M_PI / 2) < 1E-10) {
    if (fabs(point->getX() - _start.getX()) < 1E-10)
      return true;
    else
      return false;
  }

  /* If the track is not vertical */
  else {
    m = sin(_phi) / cos(_phi);
    dist = point->getY() - (_start.getY() + m * (point->getX()-_start.getX()));

    /* Use point-slope formula */
    if (fabs(dist) < 1e-10)
      return true;
    else
      return false;
  }
}


/**
 * @brief Deletes each of this Track's segments.
 */
void Track::clearSegments() {
  _segments.clear();
}


/**
 * @brief Convert this Track's attributes to a character array.
 * @details The character array returned includes the Track's starting and
 *          ending coordinates, the azimuthal angle and azimuthal weight.
 * @return a character array of this Track's attributes
 */
std::string Track::toString() {
  std::stringstream string;
  string << "Track: start, x = " << _start.getX() << ", y = " <<
    _start.getY() << " end, x = " << _end.getX() << ", y = "
                  << _end.getY() << ", phi = " << _phi;
  return string.str();
}
