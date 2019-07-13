#include "Point.h"


/**
 * @brief Constructor initializes an empty Point.
 */
Point::Point() {
  _xyz[0] = 0.0;
  _xyz[1] = 0.0;
  _xyz[2] = 0.0;
}


/**
 * @brief Destructor.
 */
Point::~Point() {
}


/**
 * @brief Converts this Point to a character representation of its attributes.
 * @details The character array includes the x-coordinate, y-coordinate, 
 *          and z-coordinate
 * @return a character array of this Point's attributes
 */
std::string Point::toString() {
  std::stringstream string;
  string << "Point: x = " << _xyz[0] << ", y = " << _xyz[1] << ", z = " << _xyz[2];
  return string.str();
}
