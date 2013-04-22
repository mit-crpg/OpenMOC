#include "Point.h"


/**
 * @brief Constructor initializes an empty point.
 */
Point::Point() { }


/**
 * @brief Destructor
 */
Point::~Point() { }


/**
 * @brief Converts this point to a character representation of its attributes.
 * @details The character array includes the x-coordinate and y-coordinates.
 * @return a character array of this point's attributes
 */
std::string Point::toString() {
    std::stringstream string;
    string << "Point: x = " << _x << ", y = " << _y;
    return string.str();
}
