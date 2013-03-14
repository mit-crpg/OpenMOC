/*
 * Point.cpp
 *
 *  Created on: Jan 18, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#include "Point.h"


/**
 * Default constructor
 */
Point::Point() { }


/**
 * Destructor
 */
Point::~Point() { }


/**
 * Converts this point to a character representation of its attributes
 * @param a character array of this point's attributes
 */
std::string Point::toString() {
	std::stringstream string;

	string << "Point: x = " << _x << ", y = " << _y;

	return string.str();
}
