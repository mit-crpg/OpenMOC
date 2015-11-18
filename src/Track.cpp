#include "Track.h"


/*
 * @brief Constructor initializes an empty Track.
 */
Track::Track() {

  /* Initialize the periodic track index to -1, indicating it has not
   * been set */
  _periodic_track_index = -1;

  /* Initialize the reflective track index to -1, indicating it has not
   * been set */
  _reflective_track_index = -1;
}



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
 * @param start_z the z-coordinate at the starting point
 * @param end_x the x-coordinate at the ending point
 * @param end_y the y-coordinate at the ending point
 * @param end_z the z-coordinate at the ending point
 * @param phi the track's azimuthal angle (\f$ \theta \in [0, \pi] \f$)
 */
 void Track::setValues(const double start_x, const double start_y,
                       const double start_z, const double end_x,
                       const double end_y, const double end_z,
                       const double phi) {
   _start.setCoords(start_x, start_y, start_z);
   _end.setCoords(end_x, end_y, end_z);
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
 * @brief Set the index of a track in a periodic cycle.
 * @details Tracks form periodic track cycles as they traverse the geometry.
 *          Tracks can be arbitrarily decomposed into periodic track cycles
 *          and this index indicates the index in a particular cycle.
 * @param index of the track in a periodic cycle
 */
void Track::setPeriodicTrackIndex(const int index) {
  _periodic_track_index = index;
}


/**
 * @brief Set the index of a track in a reflective cycle.
 * @details Tracks form reflective track cycles as they traverse the geometry.
 *          Tracks can be arbitrarily decomposed into reflective track cycles
 *          and this index indicates the index in a particular cycle.
 * @param index of the track in a reflective cycle
 */
void Track::setReflectiveTrackIndex(const int index) {
  _reflective_track_index = index;
}


/**
 * @brief Adds a segment pointer to this Track's list of segments.
 * @details This method assumes that segments are added in order of their
 *          starting location from the Track's start point.
 * @param to_add a pointer to the segment to add
 */
void Track::addSegment(segment* to_add) {
  try {
    _segments.push_back(*to_add);
  }
  catch (std::exception &e) {
    log_printf(ERROR, "Unable to add a segment to Track");
  }
}


/**
 * @brief Removes a segment from this Track's list of segments.
 * @param index the index of the segment to remove
 */
void Track::removeSegment(int index) {
  try {
    _segments.erase(_segments.begin()+index);
  }
  catch (std::exception &e) {
    log_printf(ERROR, "Unable to remove a segment from Track");
  }
}


/**
 * @brief Inserts a segment pointer into this Track's list of segments.
 * @details This method appends the new segment directly behind another
 *          segment in the Track. This is a helper method for the
 *          TrackGenerator::splitTracks(...) routine.
 * @param index the index of the segment to insert behind in the list
 * @param segment a pointer to the segment to insert
 */
void Track::insertSegment(int index, segment* segment) {
  try {
    _segments.insert(_segments.begin()+index, *segment);
  }
  catch (std::exception &e) {
    log_printf(ERROR, "Unable to insert a segment into Track");
  }
}


/**
 * @brief Sets the direction in which the flux leaving this Track along its
 *        "forward" direction is passed.
 * @details Sets whether or not to pass the outgoing flux from this Track
 *          along its "forward" direction to the "forward" direction (false)
 *          or "reverse" direction (true) of the next Track after intersection
 *          with the geometry boundary.
 * @param next_in the "forward" (false) or "reverse (true) direction
 */
void Track::setNextIn(const bool next_in) {
  _next_in = next_in;
}


/**
 * @brief Sets the direction in which the flux leaving this Track along its
 *        "reverse" direction is passed.
 * @details Sets whether or not to pass the outgoing flux from this Track
 *          along its "reverse" direction to the "forward" direction (false)
 *          or "reverse" direction (true) of the next Track after intersection
 *          with the geometry boundary.
 * @param next_out the "forward" (false) or "reverse (true) direction
 */
void Track::setNextOut(const bool next_out) {
  _next_out = next_out;
}


/**
 * @brief Sets the boundary condition for the incoming flux along the Track's
 *        "forward" direction.
 * @details The boundaryType represents vacuum (0), reflective (1), or
 *          periodic (2) boundary conditions.
 * @param bc_in boundary condition for the incoming flux in the "forward"
 *        direction
 */
void Track::setBCIn(const boundaryType bc_in) {
  _bc_in = bc_in;
}


/**
 * @brief Sets the boundary condition for the incoming flux along the Track's
 *        "reverse" direction.
 * @details The boundaryType represents vacuum (0), reflective (1), or
 *          periodic (2) boundary conditions.
 * @param bc_out boundary condition for the incoming flux in the "reverse"
 *        direction
 */
void Track::setBCOut(const boundaryType bc_out) {
  _bc_out = bc_out;
}


/**
 * @brief Sets the track going out along this Track's "forward" direction.
 * @param track_in pointer to the Track going out in the "forward" direction
 */
void Track::setTrackIn(Track* track_in) {
  _track_in = track_in;
}


/**
 * @brief Sets the track going out along this Track's "reverse" direction.
 * @param track_out pointer to the Track going out in the "reverse" direction
 */
void Track::setTrackOut(Track* track_out) {
  _track_out = track_out;
}


/**
 * @brief Returns whether to give the outgoing flux to the "forward" (false) or
 *        "reverse" (true) direction of the next Track when traveling along
 *        this Tracks's "forward" direction.
 * @return "forward" (false) "reverse" (true) direction of outgoing Track
 */
bool Track::isNextIn() const {
  return _next_in;
}


/**
 * @brief Returns whether to give the outgoing flux to the "forward" (false) or
 *        "reverse" (true) direction of the next Track when traveling along
 *        this Track's "reverse" direction.
 * @return "forward" (false) "reverse" (true) direction of outgoing Track
 */
bool Track::isNextOut() const {
  return _next_out;
}


/**
 * @brief Returns the boundary condition for the flux along the Track's
 *        "forward" direction.
 * @return vacuum (0), reflective (1), or periodic (2) reflective
 *         boundary conditions
 */
boundaryType Track::getBCIn() const {
  return _bc_in;
}


/**
 * @brief Returns the boundary condition for the flux along the Track's
 *        "reverse" direction.
 * @return vacuum (0), reflective (1), or periodic (2) reflective
 *         boundary conditions
 */
boundaryType Track::getBCOut() const {
  return _bc_out;
}


/**
 * @brief Returns a boolean to indicate whether the outgoing flux along this
 *        Track's "forward" direction should be transferred to the outgoing
 *        Track.
 * @details The bool with be false for vacuum BCs and true for all other BCs.
 * @return bool indicating whether the flux should be passed when tracking in
 *         the "forward" direction.
 */
bool Track::getTransferFluxIn() const {

  if (_bc_in == VACUUM)
    return false;
  else
    return true;
}


/**
 * @brief Returns a boolean to indicate whether the outgoing flux along this
 *        Track's "reverse" direction should be transferred to the incoming
 *        Track.
 * @details The bool with be false for vacuum BCs and true for all other BCs.
 * @return bool indicating whether the flux should be passed when tracking in
 *         the "reverse" direction.
 */
bool Track::getTransferFluxOut() const {

  if (_bc_out == VACUUM)
    return false;
  else
    return true;
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
 * @brief Get the index of a track in a periodic cycle.
 * @return index of the track in a periodic cycle
 */
int Track::getPeriodicTrackIndex() const {
  return _periodic_track_index;
}


/**
 * @brief Get the index of a track in a reflective cycle.
 * @return index of the track in a reflective cycle
 */
int Track::getReflectiveTrackIndex() const {
  return _reflective_track_index;
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
    _start.getY() << ", z = " << _start.getZ() << ", end, x = " <<
    _end.getX() << ", y = " << _end.getY() << ", z = " << _end.getZ() <<
    ", phi = " << _phi;

  return string.str();
}
