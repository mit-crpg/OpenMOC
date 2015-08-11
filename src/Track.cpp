#include "Track.h"


/*
 * @brief Constructor initializes an empty Track.
 */
Track::Track() {

  _track_refl_fwd = NULL;
  _track_refl_bwd = NULL;
  _track_prdc_fwd = NULL;
  _track_prdc_bwd = NULL;
  _refl_fwd_fwd = true;
  _refl_bwd_fwd = false;
  _prdc_fwd_fwd = true;
  _prdc_bwd_fwd = false;

  _periodic_cycle_id = -1;
  _reflective_cycle_id = -1;
  _periodic_track_index = -1;
}



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
void Track::setBCFwd(const boundaryType bc_fwd) {
  _bc_fwd = bc_fwd;
}


/**
 * @brief Sets the boundary condition for the incoming flux along the Track's
 *        "reverse" direction.
 * @details The boolean represents vacuum (false) or reflective (true)
 *          boundary conditions.
 * @param bc_out boundary condition for the incoming flux in the "reverse"
 *        direction
 */
void Track::setBCBwd(const boundaryType bc_bwd) {
  _bc_bwd = bc_bwd;
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
boundaryType Track::getBCFwd() const {
  return _bc_fwd;
}


/**
 * @brief Returns the boundary condition for the flux along the Track's
 *        "reverse" direction.
 * @return vacuum (false) or reflective (true) reflective boundary conditions
 */
boundaryType Track::getBCBwd() const {
  return _bc_bwd;
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
      log_printf(NORMAL, "Unable to add a segment to Track. Backtrace:"
                 "\n%s", e.what());
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
    log_printf(NORMAL, "Unable to remove a segment from Track");
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
    log_printf(NORMAL, "Unable to insert a segment into Track");
  }  
}


/**
 * @brief Deletes each of this Track's segments.
 */
void Track::clearSegments() {
  _segments.clear();
}


void Track::setAzimIndex(int index){
  _azim_index = index;
}


int Track::getAzimIndex() {
  return _azim_index;
}


void Track::setTrackReflFwd(Track* track){
  
  if (_refl_fwd_fwd && _end.distanceToPoint(track->getStart()) > 1.e-8)
    log_printf(NORMAL, "INCORRECT TRACK REFL FWD: (%.4f, %.4f, %.4f) -> "
               "(%.4f, %.4f, %.4f)", _end.getX(), _end.getY(), _end.getZ(),
               track->getStart()->getX(), track->getStart()->getY(),
               track->getStart()->getZ());
  if (!_refl_fwd_fwd && _end.distanceToPoint(track->getEnd()) > 1.e-8)
    log_printf(NORMAL, "INCORRECT TRACK REFL FWD: (%.4f, %.4f, %.4f) -> "
               "(%.4f, %.4f, %.4f)", _end.getX(), _end.getY(), _end.getZ(),
               track->getEnd()->getX(), track->getEnd()->getY(),
               track->getEnd()->getZ());
    
  _track_refl_fwd = track;
}


void Track::setTrackPrdcFwd(Track* track){

  if (_prdc_fwd_fwd) {
    double dx = fabs(_end.getX() - track->getStart()->getX());
    double dy = fabs(_end.getY() - track->getStart()->getY());
    double dz = fabs(_end.getZ() - track->getStart()->getZ());
    if ((dx > 1.e-8 && dy > 1.e-8) || (dx > 1.e-8 && dz > 1.e-8) ||
        (dy > 1.e-8 && dz > 1.e-8))
      log_printf(NORMAL, "INCORRECT TRACK PRDC FWD: (%.4f, %.4f, %.4f) -> "
                 "(%.4f, %.4f, %.4f)", _end.getX(), _end.getY(), _end.getZ(),
                 track->getStart()->getX(), track->getStart()->getY(),
                 track->getStart()->getZ());
  }
  else {
    double dx = fabs(_end.getX() - track->getEnd()->getX());
    double dy = fabs(_end.getY() - track->getEnd()->getY());
    double dz = fabs(_end.getZ() - track->getEnd()->getZ());
    if ((dx > 1.e-8 && dy > 1.e-8) || (dx > 1.e-8 && dz > 1.e-8) ||
        (dy > 1.e-8 && dz > 1.e-8))
      log_printf(NORMAL, "INCORRECT TRACK PRDC FWD: (%.4f, %.4f, %.4f) -> "
                 "(%.4f, %.4f, %.4f)", _end.getX(), _end.getY(), _end.getZ(),
                 track->getEnd()->getX(), track->getEnd()->getY(),
                 track->getEnd()->getZ());
  }
  
  _track_prdc_fwd = track;
}


void Track::setTrackReflBwd(Track* track){

  if (_refl_bwd_fwd && _start.distanceToPoint(track->getStart()) > 1.e-8)
    log_printf(NORMAL, "INCORRECT TRACK REFL BWD: (%.4f, %.4f, %.4f) -> "
               "(%.4f, %.4f, %.4f)", _start.getX(), _start.getY(),
               _start.getZ(), track->getStart()->getX(),
               track->getStart()->getY(), track->getStart()->getZ());
  else if (!_refl_bwd_fwd && _start.distanceToPoint(track->getEnd()) > 1.e-8)
    log_printf(NORMAL, "INCORRECT TRACK REFL BWD: (%.4f, %.4f, %.4f) -> "
               "(%.4f, %.4f, %.4f)", _start.getX(), _start.getY(),
               _start.getZ(), track->getEnd()->getX(),
               track->getEnd()->getY(), track->getEnd()->getZ());
  
  _track_refl_bwd = track;
}


void Track::setTrackPrdcBwd(Track* track){

  if (_prdc_bwd_fwd) {
    double dx = fabs(_start.getX() - track->getStart()->getX());
    double dy = fabs(_start.getY() - track->getStart()->getY());
    double dz = fabs(_start.getZ() - track->getStart()->getZ());
    if ((dx > 1.e-8 && dy > 1.e-8) || (dx > 1.e-8 && dz > 1.e-8) ||
        (dy > 1.e-8 && dz > 1.e-8))
      log_printf(NORMAL, "INCORRECT TRACK PRDC BWD: (%.4f, %.4f, %.4f) -> "
                 "(%.4f, %.4f, %.4f)", _start.getX(), _start.getY(),
                 _start.getZ(), track->getStart()->getX(),
                 track->getStart()->getY(), track->getStart()->getZ());
  }
  else {
    double dx = fabs(_start.getX() - track->getEnd()->getX());
    double dy = fabs(_start.getY() - track->getEnd()->getY());
    double dz = fabs(_start.getZ() - track->getEnd()->getZ());
    if ((dx > 1.e-8 && dy > 1.e-8) || (dx > 1.e-8 && dz > 1.e-8) ||
        (dy > 1.e-8 && dz > 1.e-8))
      log_printf(NORMAL, "INCORRECT TRACK PRDC BWD: (%.4f, %.4f, %.4f) -> "
                 "(%.4f, %.4f, %.4f)", _start.getX(), _start.getY(),
                 _start.getZ(), track->getEnd()->getX(),
                 track->getEnd()->getY(), track->getEnd()->getZ());
  }

  _track_prdc_bwd = track;
}


Track* Track::getTrackReflFwd() {
  return _track_refl_fwd;
}


Track* Track::getTrackReflBwd() {
  return _track_refl_bwd;
}


Track* Track::getTrackPrdcFwd() {
  return _track_prdc_fwd;
}


Track* Track::getTrackPrdcBwd() {
  return _track_prdc_bwd;
}


void Track::setXYIndex(int index) {
  _xy_index = index;
}


int Track::getXYIndex() {
  return _xy_index;
}


void Track::setReflFwdFwd(bool fwd){
  _refl_fwd_fwd = fwd;
}


void Track::setReflBwdFwd(bool fwd){
  _refl_bwd_fwd = fwd;
}


void Track::setPrdcFwdFwd(bool fwd){
  _prdc_fwd_fwd = fwd;
}


void Track::setPrdcBwdFwd(bool fwd){
  _prdc_bwd_fwd = fwd;
}


bool Track::getReflFwdFwd(){
  return _refl_fwd_fwd;
}


bool Track::getReflBwdFwd(){
  return _refl_bwd_fwd;
}


bool Track::getPrdcFwdFwd(){
  return _prdc_fwd_fwd;
}


bool Track::getPrdcBwdFwd(){
  return _prdc_bwd_fwd;
}


void Track::setPeriodicCycleId(int id) {
  _periodic_cycle_id = id;
}


void Track::setPeriodicTrackIndex(int index) {
  _periodic_track_index = index;
}

 
int Track::getPeriodicCycleId() {
  return _periodic_cycle_id;
}


int Track::getPeriodicTrackIndex() {
  return _periodic_track_index;
}


void Track::setReflectiveCycleId(int id) {
  _reflective_cycle_id = id;
}

 
int Track::getReflectiveCycleId() {
  return _reflective_cycle_id;
}
