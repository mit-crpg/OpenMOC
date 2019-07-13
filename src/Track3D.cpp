#include "Track3D.h"


/**
 * @brief Constructor initializes Track3D with negative values and indexes.
 */
Track3D::Track3D() : Track() {

  _theta = -1.;
  _polar_index = -1;
  _z_index = -1;
  _lz_index = -1;
  _cycle_index = -1;
  _cycle_track_index = -1;
  _train_index = -1;
  _cycle_fwd = false;
}


/**
 * @brief Destructor.
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
 * @param theta The polar angle
 */
void Track3D::setTheta(const double theta) {
  _theta = theta;
}


/**
 * @brief Set the index for the Track's polar angle index.
 * @details The polar angle index corresponds to an array of all
 *          polar angles for \f$ \theta \in [0, \pi] \f$ owned by
 *          the TrackGenerator class.
 * @param index The polar angle index
 */
void Track3D::setPolarIndex(int index) {
  _polar_index = index;
}


/**
 * @brief Set the Track's z index into the Track3D stacks array.
 * @param index The Track's z index.
 */
void Track3D::setZIndex(int index) {
  _z_index = index;
}


/**
 * @brief Set the Track's lz index into the Track3D cycles array.
 * @param index The Track's lz index.
 */
void Track3D::setLZIndex(int index) {
  _lz_index = index;
}


/**
 * @brief Set the Track's cycle index into the Track3D cycles array.
 * @param index The Track's cycle index.
 */
void Track3D::setCycleIndex(int index) {
  _cycle_index = index;
}


/**
 * @brief Set the index into the reflective Track cycle that this
 *        Track lies above.
 * @param index The Track's index into the reflective Track cycle that this
 *        Track lies above.
 */
void Track3D::setCycleTrackIndex(int index) {
  _cycle_track_index = index;
}


/**
 * @brief Set the Track's last index into the Track3D cycles array.
 * @param index The Track's last index into the Track3D cycles array.
 */
void Track3D::setTrainIndex(int index) {
  _train_index = index;
}



/**
 * @brief Return the Track's polar angle (with respect to the positive z-axis).
 * @return the polar angle \f$ \theta \in [0, \frac{\pi}{2}] \f$
 */
double Track3D::getTheta() const {
  return _theta;
}


/**
 * @brief Get the index for the Track's polar angle index.
 * @details The polar angle index corresponds to an array of all
 *          polar angles for \f$ \theta \in [0, \pi] \f$ owned by
 *          the TrackGenerator class.
 * @return The polar angle index
 */
int Track3D::getPolarIndex() {
  return _polar_index;
}


/**
 * @brief Get the Track's z index into the Track3D stacks array.
 * @return The Track's z index.
 */
int Track3D::getZIndex() {
  return _z_index;
}


/**
 * @brief Get the Track's lz index into the Track3D cycles array.
 * @return The Track's ;z index.
 */
int Track3D::getLZIndex() {
  return _lz_index;
}


/**
 * @brief Get the Track's cycle index into the Track3D cycles array.
 * @return The Track's cycle index.
 */
int Track3D::getCycleIndex() {
  return _cycle_index;
}


/**
 * @brief Get the index into the reflective Track cycle that this
 *        Track lies above.
 * @return The Track's index into the reflective Track cycle that this
 *         Track lies above.
 */
int Track3D::getCycleTrackIndex() {
  return _cycle_track_index;
}


/**
 * @brief Get the Track's last index into the Track3D cycles array.
 * @return The Track's last index into the Track3D cycles array.
 */
int Track3D::getTrainIndex() {
  return _train_index;
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


/**
 * @brief Set the values for the Track's start and end point.
 * @param x0 the x-coordinate at the starting point
 * @param y0 the y-coordinate at the starting point
 * @param z0 the z-coordinate at the starting point
 * @param x1 the x-coordinate at the ending point
 * @param y1 the y-coordinate at the ending point
 * @param z1 the z-coordinate at the ending point
 */
void Track3D::setCoords(double x0, double y0, double z0,
                        double x1, double y1, double z1) {
  _start.setCoords(x0, y0, z0);
  _end.setCoords(x1, y1, z1);
}


/**
 * @brief Set a boolean indicating whether the track cycle is going forward
 *        along the Track's forward direction.
 * @param fwd A boolean indicating whether the track cycle is going forward
 *        along the Track's forward direction.
 */
void Track3D::setCycleFwd(bool fwd) {
  _cycle_fwd = fwd;
}


/**
 * @brief Get a boolean indicating whether the track cycle is going forward
 *        along the Track's forward direction.
 * @return A boolean indicating whether the track cycle is going forward
 *         along the Track's forward direction.
 */
bool Track3D::getCycleFwd() {
  return _cycle_fwd;
}
