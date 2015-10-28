#include "Track2D.h"


/*
 * @brief Constructor initializes an empty Track2D.
 */
Track2D::Track2D() : Track() { }



/**
 * @brief Destructor clears the Track segments container.
 */
Track2D::~Track2D() {
}


/**
 * @brief Set the values for the Track's start and end point and angle.
 * @param start_x the x-coordinate at the starting point
 * @param start_y the y-coordinate at the starting point
 * @param end_x the x-coordinate at the ending point
 * @param end_y the y-coordinate at the ending point
 * @param phi the track's azimuthal angle (\f$ \phi \in [0, 2 \pi] \f$)
 */
void Track2D::setValues(const double start_x, const double start_y,
                        const double end_x, const double end_y,
                        const double phi) {
   _start.setCoords(start_x, start_y);
   _end.setCoords(end_x, end_y);
   _phi = phi;
}


/**
 * @brief Convert this Track's attributes to a character array.
 * @details The character array returned includes the Track's starting and
 *          ending coordinates, the azimuthal angle and azimuthal weight.
 * @return a character array of this Track's attributes
 */
std::string Track2D::toString() {
  std::stringstream string;
  string << "Track2D: start, x = " << _start.getX() << ", y = " <<
    _start.getY() << ", z = " << ", end, x = " <<
    _end.getX() << ", y = " << _end.getY() <<
    ", phi = " << _phi;
  return string.str();
}


/**
 * @brief Set the values for the Track's start and end point.
 * @param x0 the x-coordinate at the starting point
 * @param y0 the y-coordinate at the starting point
 * @param x1 the x-coordinate at the ending point
 * @param y1 the y-coordinate at the ending point
 */
void Track2D::setCoords(double x0, double y0,
                        double x1, double y1){
  _start.setCoords(x0, y0);
  _end.setCoords(x1, y1);
}
