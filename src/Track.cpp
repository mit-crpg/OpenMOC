#include "Track.h"


/**
 * @brief Constructor initializes an empty Track.
 */
Track::Track() {

  /* Initialize the pointers to reflective and periodic tracks to NULL */
  _track_next_fwd = -1;
  _track_next_bwd = -1;
  _track_prdc_fwd = -1;
  _track_prdc_bwd = -1;
  _track_refl_fwd = -1;
  _track_refl_bwd = -1;

  /* Initialize booleans indicating whether the reflective tracks in the
   * forward and backward direction point enter the track in the forward
   * direction */
  _next_fwd_fwd = true;
  _next_bwd_fwd = false;

  /* Initialize the cycle ids and periodic track index to -1 (not set) */
  _num_segments = 0;
  _surface_in = -1;
  _surface_out = -1;
  _domain_fwd = -1;
  _domain_bwd = -1;
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
 * @param end_x the x-coordinate at the ending point
 * @param end_y the y-coordinate at the ending point
 * @param phi the track's azimuthal angle (\f$ \phi \in [0, 2 \pi] \f$)
 */
void Track::setValues(const double start_x, const double start_y,
                      const double end_x, const double end_y,
                      const double phi) {
   _start.setCoords(start_x, start_y);
   _end.setCoords(end_x, end_y);
   _phi = phi;
}


/**
 * @brief Set the values for the Track's start and end point.
 * @param x0 the x-coordinate at the starting point
 * @param y0 the y-coordinate at the starting point
 * @param x1 the x-coordinate at the ending point
 * @param y1 the y-coordinate at the ending point
 */
void Track::setCoords(double x0, double y0,
                      double x1, double y1) {
  _start.setCoords(x0, y0);
  _end.setCoords(x1, y1);
}


/**
 * @brief Initializes a Track's unique ID.
 * @details This is set by the trackgenerator to correspond to the Track's
 *          location in a 2D ragged array of all tracks.
 * @param uid The Track's unique ID
 */
void Track::setUid(int uid) {
  _uid = uid;
}

/**
 * @brief Set the Track's azimuthal angle.
 * @param phi The azimuthal angle
 */
void Track::setPhi(const double phi) {
  _phi = phi;
}


/**
 * @brief Sets the boundary condition for the incoming flux along the Track's
 *        "forward" direction.
 * @details The boolean represents vacuum (false) or reflective (true)
 *          boundary conditions.
 * @param bc_fwd Boundary condition for the incoming flux in the "forward"
 *        direction
 */
void Track::setBCFwd(const boundaryType bc_fwd) {
  _bc_fwd = bc_fwd;
}


/**
 * @brief Sets the boundary condition for the incoming flux along the Track's
 *        "reverse" direction.
 * @details The boolean represents vacuum (false) or reflective (true)
 *          boundary conditions.
 * @param bc_bwd Boundary condition for the incoming flux in the "reverse"
 *        direction
 */
void Track::setBCBwd(const boundaryType bc_bwd) {
  _bc_bwd = bc_bwd;
}


/**
 * @brief Set the surface at which the Track starts.
 * @param surface_in surface from which the Track originates in the domain.
 */
void Track::setSurfaceIn(int surface_in) {
  _surface_in = surface_in;
}


/**
 * @brief Set the surface at which the Track ends.
 * @param surface_out surface at which the Track ends in the domain.
 */
void Track::setSurfaceOut(int surface_out) {
  _surface_out = surface_out;
}


/**
 * @brief Set the domain at the end of the Track in the forward direction.
 * @param neighbor Neighbor domain at the end of the Track
 */
void Track::setDomainFwd(int neighbor) {
  _domain_fwd = neighbor;
}


/**
 * @brief Set the domain at the end of the Track in the backward direction.
 * @param neighbor Neighbor domain at the beginning of the Track
 */
void Track::setDomainBwd(int neighbor) {
  _domain_bwd = neighbor;
}


/**
 * @brief Return the Track's length.
 * @return The Track's length
 */
double Track::getLength() {
  return _start.distanceToPoint(&_end);
}


/**
 * @brief Get the surface from which the track originates in the domain.
 * @return The id of the Surface at the beginning of the Track
 */
int Track::getSurfaceIn() {
    return _surface_in;
}


/**
 * @brief Get the surface at which the Track ends in the domain.
 * @return The id of the Surface at the end of the Track
 */
int Track::getSurfaceOut() {
    return _surface_out;
}


/**
 * @brief Get the domain connected at the front end of the Track.
 * @return The id of the domain in front of the Track
 */
int Track::getDomainFwd() {
  return _domain_fwd;
}


/**
 * @brief Get the domain connected at the back end of the Track.
 * @return The id of the domain behind the Track
 */
int Track::getDomainBwd() {
  return _domain_bwd;
}


/**
 * @brief Adds a segment to this Track's list of segments.
 * @details This method assumes that segments are added in order of their
 *          starting location from the Track's start point.
 * @param segment A pointer to the segment
 */
void Track::addSegment(segment* segment) {

  try {
    _segments.push_back(*segment);
  }
  catch (std::exception &e) {
      log_printf(NORMAL, "Unable to add a segment to Track. Backtrace:"
                 "\n%s", e.what());
  }
}


/**
 * @brief Removes a segment from this Track's list of segments.
 * @param index The index of the segment to remove
 */
void Track::removeSegment(int index) {
  try {
    _segments.erase(_segments.begin()+index);
  }
  catch (std::exception &e) {
    log_printf(NORMAL, "Unable to remove a segment from Track");
  }
}


/**
 * @brief Inserts a segment struct into this Track's list of segments.
 * @details This method appends the new segment directly behind another
 *          segment in the Track. This is a helper method for the
 *          TrackGenerator::splitTracks(...) routine.
 * @param index The index of the segment to insert behind in the list
 * @param segment A pointer to the segment to insert
 */
void Track::insertSegment(int index, segment* segment) {
  try {
    _segments.insert(_segments.begin()+index, *segment);
  }
  catch (std::exception &e) {
    log_printf(NORMAL, "Unable to insert a segment into Track");
  }
}


/**
 * @brief Deletes each of this Track's segments.
 */
void Track::clearSegments() {
  _segments.clear();
}


/**
 * @brief Set a Track's azimuthal angle index.
 * @param index The azimuthal angle index
 */
void Track::setAzimIndex(int index) {
  _azim_index = index;
}


/**
 * @brief Set a Track's link index.
 * @details In generating 3D tracks, we need to know the tracks
 *          linking this track with both periodic and reflective
 *          boundary conditions. Therefore, we generate a periodic
 *          chain of tracks that extend periodically from the y-min
 *          to the y-max boundary. This is the index of each track in
 *          its' periodic chain.
 * @param index The link index
 */
void Track::setLinkIndex(int index) {
  _link_index = index;
}


/**
 * @brief Get a Track's azimuthal angle index.
 * @return The azimuthal angle index
 */
int Track::getAzimIndex() {
  return _azim_index;
}


/**
 * @brief Get a Track's link index.
 * @details In generating 3D tracks, we need to know the tracks
 *          linking this track with both periodic and reflective
 *          boundary conditions. Therefore, we generate a periodic
 *          chain of tracks that extend periodically from the y-min
 *          to the y-max boundary. This is the index of each track in
 *          its' periodic chain.
 * @return The link index
 */
int Track::getLinkIndex() {
  return _link_index;
}


/**
 * @brief Set a pointer to the reflective Track in the forward direction.
 * @param track_id the id of the reflective track in the forward direction
 */
void Track::setTrackNextFwd(long track_id) {
  _track_next_fwd = track_id;
}


/**
 * @brief Set a pointer to the reflective Track in the backward direction.
 * @param track_id the id of the reflective track in the backward direction
 */
void Track::setTrackNextBwd(long track_id) {
  _track_next_bwd = track_id;
}


/**
 * @brief Set a pointer to the periodic Track in the forward direction.
 * @param track_id the id of the periodic track in the forward direction
 */
void Track::setTrackPrdcFwd(long track_id) {
  _track_prdc_fwd = track_id;
}


/**
 * @brief Set a pointer to the periodic Track in the backward direction.
 * @param track_id the id of the periodic track in the backward direction
 */
void Track::setTrackPrdcBwd(long track_id) {
  _track_prdc_bwd = track_id;
}


/**
 * @brief Set a pointer to the reflective Track in the forward direction.
 * @param track_id the id of the reflective track in the forward direction
 */
void Track::setTrackReflFwd(long track_id) {
  _track_refl_fwd = track_id;
}


/**
 * @brief Set a pointer to the reflective Track in the backward direction.
 * @param track_id the id of the reflective track in the backward direction
 */
void Track::setTrackReflBwd(long track_id) {
  _track_refl_bwd = track_id;
}


/**
 * @brief Get a pointer to the reflective Track in the forward direction.
 * @return the id of the reflective track in the forward direction
 */
long Track::getTrackNextFwd() {
  return _track_next_fwd;
}


/**
 * @brief Get a pointer to the reflective Track in the backward direction.
 * @return the id of the reflective track in the backward direction
 */
long Track::getTrackNextBwd() {
  return _track_next_bwd;
}


/**
 * @brief Get a pointer to the periodic Track in the forward direction.
 * @return the id of the periodic track in the forward direction
 */
long Track::getTrackPrdcFwd() {
  return _track_prdc_fwd;
}


/**
 * @brief Get a pointer to the periodic Track in the backward direction.
 * @return the id of the periodic track in the backward direction
 */
long Track::getTrackPrdcBwd() {
  return _track_prdc_bwd;
}


/**
 * @brief Get a pointer to the reflective Track in the forward direction.
 * @return the id of the reflective track in the forward direction
 */
long Track::getTrackReflFwd() {
  return _track_refl_fwd;
}


/**
 * @brief Get a pointer to the reflective Track in the backward direction.
 * @return the id of the reflective track in the backward direction
 */
long Track::getTrackReflBwd() {
  return _track_refl_bwd;
}


/**
 * @brief Set the xy index of this Track.
 * @param index The xy index of this Track
 */
void Track::setXYIndex(int index) {
  _xy_index = index;
}


/**
 * @brief Get the xy index of this Track.
 * @return The xy index of this Track
 */
int Track::getXYIndex() {
  return _xy_index;
}


/**
 * @brief Set whether the reflective track in the forward direction is pointing
 *        in forward direction.
 * @param fwd Boolean indicating whether reflective track in the forward
 *        direction is point in forward direction
 */
void Track::setNextFwdFwd(bool fwd) {
  _next_fwd_fwd = fwd;
}


/**
 * @brief Set whether the reflective track in the backward direction is pointing
 *        in forward direction.
 * @param fwd Boolean indicating whether reflective track in the backward
 *        direction is point in forward direction
 */
void Track::setNextBwdFwd(bool fwd) {
  _next_bwd_fwd = fwd;
}


/**
 * @brief Get whether the reflective track in the forward direction is pointing
 *        in forward direction.
 * @return Boolean indicating whether reflective track in the forward
 *         direction is point in forward direction
 */
bool Track::getNextFwdFwd() {
  return _next_fwd_fwd;
}


/**
 * @brief Get whether the reflective track in the backward direction is pointing
 *        in forward direction.
 * @return Boolean indicating whether reflective track in the backward
 *         direction is point in forward direction
 */
bool Track::getNextBwdFwd() {
  return _next_bwd_fwd;
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
    _start.getY() << ", end, x = " <<
    _end.getX() << ", y = " << _end.getY() <<
    ", phi = " << _phi;
  return string.str();
}
