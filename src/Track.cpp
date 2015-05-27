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
 * @return the azimuthal angle \f$ \phi \in [0, \pi] \f$
 */
double Track::getPhi() const {
  return _phi;
}


double Track::getLength() {
  return _start.distanceToPoint(&_end);
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
 * @brief Deletes each of this Track's segments.
 */
void Track::clearSegments() {
  _segments.clear();
}


void Track::setAzimIndex(int azim_index){
  _azim_index = azim_index;
}


int Track::getAzimIndex(){
  return _azim_index;
}
