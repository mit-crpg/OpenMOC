#include "Track3D.h"


/*
 * @brief Constructor initializes an empty Track3D.
 */
Track3D::Track3D() : Track() { }



/**
 * @brief Destructor clears the Track segments container.
 */
Track3D::~Track3D() {
}


/**
 * @brief Set the values for the Track's start and end point and angle.
 * @param start_x the x-coordinate at the starting point
 * @param start_y the y-coordinate at the starting point
 * @param start_z the z-coordinate at the starting point
 * @param end_x the x-coordinate at the ending point
 * @param end_y the y-coordinate at the ending point
 * @param end_z the z-coordinate at the ending point
 * @param phi the track's azimuthal angle (\f$ \phi \in [0, 2 \pi] \f$)
 * @param theta the track's polar angle (\f$ \theta \in [0, \pi] \f$)
 */
void Track3D::setValues(const double start_x, const double start_y,
                        const double start_z, const double end_x,
                        const double end_y, const double end_z,
                        const double phi, const double theta) {
   _start.setCoords(start_x, start_y, start_z);
   _end.setCoords(end_x, end_y, end_z);
   _phi = phi;
   _theta = theta;
}


/**
 * @brief Set the Track's polar angle.
 * @param theta the polar angle
 */
void Track3D::setTheta(const double theta) {
  _theta = theta;
}


/**
 * @brief Set the index for the Track's polar angle index.
 * @details The polar angle index corresponds to a an array of all
 *          polar angles for \f$ \theta \in [0, \pi] \f$ owned by
 *          the TrackGenerator class.
 * @param index the polar angle index
 */
void Track3D::setPolarIndex(const int index) {
  _polar_index = index;
}


void Track3D::setLZIndex(const int index) {
  _lz_index = index;
}


void Track3D::setTrackInPolarIndex(const int index) {
  _track_in_polar_index = index;
}


void Track3D::setTrackInLZIndex(const int index) {
  _track_in_lz_index = index;
}



void Track3D::setTrackOutPolarIndex(const int index) {
  _track_out_polar_index = index;
}


void Track3D::setTrackOutLZIndex(const int index) {
  _track_out_lz_index = index;
}


/**
 * @brief Return the Track's polar angle (with respect to the positive z-axis).
 * @return the polar angle \f$ \theta \in [0, \frac{\pi}{2}] \f$
 */
double Track3D::getTheta() const {
  return _theta;
}


int Track3D::getPolarIndex() const {
  return _polar_index;
}


int Track3D::getLZIndex() const {
  return _lz_index;
}


int Track3D::getTrackInPolarIndex() const {
  return _track_in_polar_index;
}


int Track3D::getTrackInLZIndex() const {
  return _track_in_lz_index;
}


int Track3D::getTrackOutPolarIndex() const {
  return _track_out_polar_index;
}


int Track3D::getTrackOutLZIndex() const {
  return _track_out_lz_index;
}


/**
 * @brief Convert this Track's attributes to a character array.
 * @details The character array returned includes the Track's starting and
 *          ending coordinates, the azimuthal angle and azimuthal weight.
 * @return a character array of this Track's attributes
 */
std::string Track3D::toString() {
  std::stringstream string;
  string << "Track3D: start, x = " << _start.getX() << ", y = " <<
    _start.getY() << ", z = " << _start.getZ() << " end, x = " <<
    _end.getX() << ", y = " << _end.getY() << ", z = " << _end.getZ() <<
    ", phi = " << _phi << ", theta = " << _theta;
  return string.str();
}


void Track3D::setOutgoingTrack(int polar, int lz, int train){
  _track_out_polar_index = polar;
  _track_out_lz_index = lz;
  _track_out_train_index = train;
}


void Track3D::setIncomingTrack(int polar, int lz, int train){
  _track_in_polar_index = polar;
  _track_in_lz_index = lz;
  _track_in_train_index = train;
}


void Track3D::setCoords(double x0, double y0, double z0,
                        double x1, double y1, double z1){
  _start.setCoords(x0, y0, z0);
  _end.setCoords(x1, y1, z1);
}
