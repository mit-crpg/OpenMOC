/*
 * Surface.cpp
 *
 *  Created on: Jan 9, 2012
 *      Author: wbinventor
 */

#include "Surface.h"

/* _n keeps track of the number of surfaces instantiated */
short int Surface::_n = 0;


/**
 * Default Surface constructor
 * @param id the surface id
 * @param type the surface type
 * @param boundary this surface's boundary type
 */
Surface::Surface(const short int id, const surfaceType type,
					const boundaryType boundary){
	_uid = _n;
	_id = id;
	_type = type;
	_n++;
	_boundary = boundary;
}


/**
 * Surface Destructor
 */
Surface::~Surface() { }


/**
 * Return the surface's uid
 * @return the surface's uid
 */
short int Surface::getUid() const {
	return _uid;
}


/**
 * Return the surface's id
 * @return the surface's id
 */
short int Surface::getId() const {
    return _id;
}


/**
 * Return the surface's type
 * @return the surface type
 */
surfaceType Surface::getType() const {
    return _type;
}


/**
 * Returns the surface's boundary type
 * @return the boundary type (REFLECTIVE, VACUUM or BOUNDARY_NONE)
 */
boundaryType Surface::getBoundary(){
       return _boundary;
}


/**
 * Return true or false if a point is on or off of a surface.
 * @param point pointer to the point of interest
 * @return true (on) or false (off)
 */
bool Surface::onSurface(Point* point) {

	/* Uses a threshold to determine whether the point is on the surface */
	if (abs(evaluate(point)) < ON_SURFACE_THRESH)
		return true;
	else
		return false;
}


/**
 * Return true or false if a localcoord is on or off of a surface
 * @param point pointer to the localcoord of interest
 * @return true (on) or false (off)
 */
bool Surface::onSurface(LocalCoords* coord) {
	return onSurface(coord->getPoint());
}


/**
 * Plane constructor
 * @param id the surface id
 * @param boundary this surface's boundary type
 * @param A the first coefficient in A * x + B * y + C = 0;
 * @param B the second coefficient in A * x + B * y + C = 0;
 * @param C the third coefficient in A * x + B * y + C = 0;
 */
Plane::Plane(const short int id, const boundaryType boundary,
	     const double A, const double B,
	     const double C, const double D): Surface(id, PLANE, boundary) {
	_A = A;
	_B = B;
	_C = C;
	_D = D;
}


/**
 * Converts this Plane's attributes to a character array
 * @return a character array of this plane's attributes
 */
std::string Plane::toString() {
	std::stringstream string;

	string << "Surface id = " << _id << ", type = PLANE " << ", A = "
			<< _A << ", B = " << _B << ", C = " << _C << ", D = " << _D;

	return string.str();
}


/**
* Finds the intersection point with this plane from a given point and
* trajectory defined by an angle
* @param point pointer to the point of interest
* @param angle the angle defining the trajectory in radians
* @param points pointer to a point to store the intersection point
* @return the number of intersection points (0 or 1)
*/
inline int Plane::intersection(Point* point, double angle, Point* points) {

	double x0 = point->getX();
	double y0 = point->getY();

	int num = 0; 			/* number of intersections */
	double xcurr, ycurr;	/* coordinates of current intersection point */

	/* The track is vertical */
	if ((fabs(angle - (M_PI / 2))) < 1.0e-10) {

		/* The plane is also vertical => no intersections */
		if (_B == 0)
			return 0;

		/* The plane is not vertical */
		else {
			xcurr = x0;
			ycurr = (-_A * x0 - _D) / _B;
			points->setCoords(xcurr, ycurr);
			/* Check that point is in same direction as angle */
			if (angle < M_PI && ycurr > y0)
				num++;
			else if (angle > M_PI && ycurr < y0)
				num++;
			return num;
		}
	}

	/* If the track isn't vertical */
	else {
		double m = sin(angle) / cos(angle);

		/* The plane and track are parallel, no intersections */
		if (fabs(-_A/_B - m) < 1e-11 && _B != 0)
			return 0;

		else {
			xcurr = -(_B * (y0 - m * x0) + _D)
					/ (_A + _B * m);
			ycurr = y0 + m * (xcurr - x0);
			points->setCoords(xcurr, ycurr);
			if (angle < M_PI && ycurr > y0)
				num++;
			else if (angle > M_PI && ycurr < y0)
				num++;

			return num;
		}
	}
}


/**
 * Returns the minimum x value on this surface
 * @return the minimum x value
 */
double Plane::getXMin(){
	log_printf(ERROR, "Plane::getXMin not implemented");
	return std::numeric_limits<double>::infinity();
}


/**
 * Returns the maximum x value on this surface
 * @return the maximum x value
 */
double Plane::getXMax(){
	log_printf(ERROR, "Plane::getXMax not implemented");
	return std::numeric_limits<double>::infinity();
}


/**
 * Returns the minimum y value on this surface
 * @return the minimum y value
 */
double Plane::getYMin(){
	log_printf(ERROR, "Plane::getYMin not implemented");
	return std::numeric_limits<double>::infinity();
}


/**
 * Returns the maximum y value on this surface
 * @return the maximum y value
 */
double Plane::getYMax(){
	log_printf(ERROR, "Plane::getYMax not implemented");
	return std::numeric_limits<double>::infinity();
}


/**
 * XPlane constructor for a plane parallel to the x-axis
 * @param id the surface id
 * @param boundary this surface's boundary type
 * @param C the location of the plane along the y-axis
 */
XPlane::XPlane(const short int id, const boundaryType boundary,
		const double D): Plane(id, boundary, 0, 1, 0, -D) {
	_type = XPLANE;
}


/**
 * Converts this XPlane's attributes to a character array
 * @param a character array of this plane's attributes
 */
std::string XPlane::toString() {
	std::stringstream string;

	string << "Surface id = " << _id << ", type = XPLANE " << ", A = "
			<< _A << ", B = " << _B << ", C = " << _C << ", D = " << _D;

	return string.str();
}


/**
 * Returns the minimum x value on this surface
 * @return the minimum x value
 */
double XPlane::getXMin(){
	return -_D;
}


/**
 * Returns the maximum x value on this surface
 * @return the maximum x value
 */
double XPlane::getXMax(){
	return -_D;
}


/**
 * Returns the minimum y value on this surface
 * @return the minimum y value
 */
double XPlane::getYMin(){
	return std::numeric_limits<double>::infinity();
}


/**
 * Returns the maximum y value on this surface
 * @return the maximum y value
 */
double XPlane::getYMax(){
	return std::numeric_limits<double>::infinity();
}


/**
 * YPlane constructor for a plane parallel to the y-axis
 * @param id the surface id
 * @param boundary the surface boundary type
 * @param D the location of the plane along the x-axis
 */
YPlane::YPlane(const short int id, const boundaryType boundary,
		const double D): Plane(id, boundary, 1, 0, 0, -D) {
	_type = YPLANE;
}


/**
 * Converts this YPlane's attributes to a character array
 * @param a character array of this plane's attributes
 */
std::string YPlane::toString() {
	std::stringstream string;

	string << "Surface id = " << _id << ", type = YPLANE " << ", A = "
			<< _A << ", B = " << _B << ", C = " << _C << ", D = " << _D;

	return string.str();
}


/**
 * Returns the minimum x value on this surface
 * @return the minimum x value
 */
double YPlane::getXMin(){
	return std::numeric_limits<double>::infinity();

}


/**
 * Returns the maximum x value on this surface
 * @return the maximum x value
 */
double YPlane::getXMax(){
	return std::numeric_limits<double>::infinity();
}


/**
 * Returns the minimum y value on this surface
 * @return the minimum y value
 */
double YPlane::getYMin(){
	return -_D;
}


/**
 * Returns the maximum y value on this surface
 * @return the maximum y value
 */
double YPlane::getYMax(){
	return -_D;
}


/**
 * ZPlane constructor for a plane parallel to the z-axis
 * @param id the surface id
 * @param boundary the surface boundary type
 * @param D the location of the plane along the x-axis
 */
ZPlane::ZPlane(const short int id, const boundaryType boundary,
		const double D): Plane(id, boundary, 0, 0, 1, -D) {
	_type = ZPLANE;
}


/**
 * Converts this ZPlane's attributes to a character array
 * @param a character array of this plane's attributes
 */
std::string ZPlane::toString() {
	std::stringstream string;

	string << "Surface id = " << _id << ", type = ZPLANE " << ", A = "
			<< _A << ", B = " << _B << ", C = " << _C << ", D = " << _D;

	return string.str();
}


/**
 * Returns the minimum x value on this surface
 * @return the minimum x value
 */
double ZPlane::getXMin(){
	return std::numeric_limits<double>::infinity();

}


/**
 * Returns the maximum x value on this surface
 * @return the maximum x value
 */
double ZPlane::getXMax(){
	return std::numeric_limits<double>::infinity();
}


/**
 * Returns the minimum y value on this surface
 * @return the minimum y value
 */
double ZPlane::getYMin(){
	return std::numeric_limits<double>::infinity();
}


/**
 * Returns the maximum y value on this surface
 * @return the maximum y value
 */
double ZPlane::getYMax(){
	return std::numeric_limits<double>::infinity();
}


/**
 * Circle constructor
 * @param id the surface id
 * @param boundary the surface boundary type
 * @param x the x-coordinte of the circle center
 * @param y the y-coordinate of the circle center
 * @param radius the radius of the circle
 */
Circle::Circle(const short int id, const boundaryType boundary, const double x,
		const double y, const double radius): Surface(id, CIRCLE, boundary) {
	_A = 1;
	_B = 1;
	_C = -2*x;
	_D = -2*y;
	_E = x*x + y*y - radius*radius;
	_radius = radius;
	center.setX(x);
	center.setY(y);
}


/**
 * Finds the intersection point with this circle from a given point and
 * trajectory defined by an angle (0, 1, or 2 points)
 * @param point pointer to the point of interest
 * @param angle the angle defining the trajectory in radians
 * @param points pointer to a an array of points to store intersection points
 * @return the number of intersection points (0 or 1)
 */
int Circle::intersection(Point* point, double angle, Point* points) {

	double x0 = point->getX();
	double y0 = point->getY();
	double xcurr, ycurr;
	int num = 0;			/* Number of intersection points */
	double a, b, c, q, discr;

	/* If the track is vertical */
	if ((fabs(angle - (M_PI / 2))) < 1.0e-10) {
		/* Solve for where the line x = x0 and the surface F(x,y) intersect
		 * Find the y where F(x0, y) = 0
		 * Substitute x0 into F(x,y) and rearrange to put in
		 * the form of the quadratic formula: ay^2 + by + c = 0
		 */
		a = _B * _B;
		b = _D;
		c = _A * x0 * x0 + _C * x0 + _E;

		discr = b*b - 4*a*c;

		/* There are no intersections */
		if (discr < 0)
			return 0;

		/* There is one intersection (ie on the surface) */
		else if (discr == 0) {
			xcurr = x0;
			ycurr = -b / (2*a);
			points[num].setCoords(xcurr, ycurr);
			if (angle < M_PI && ycurr > y0)
				num++;
			else if (angle > M_PI && ycurr < y0)
				num++;
			return num;
		}

		/* There are two intersections */
		else {
			xcurr = x0;
			ycurr = (-b + sqrt(discr)) / (2 * a);
			points[num].setCoords(xcurr, ycurr);
			if (angle < M_PI && ycurr > y0)
				num++;
			else if (angle > M_PI && ycurr < y0)
				num++;

			xcurr = x0;
			ycurr = (-b - sqrt(discr)) / (2 * a);
			points[num].setCoords(xcurr, ycurr);
			if (angle < M_PI && ycurr > y0)
				num++;
			else if (angle > M_PI && ycurr < y0)
				num++;
			return num;
		}
	}

	/* If the track isn't vertical */
	else {
		/*Solve for where the line y-y0 = m*(x-x0) and the surface F(x,y) intersect
		 * Find the (x,y) where F(x, y0 + m*(x-x0)) = 0
		 * Substitute the point-slope formula for y into F(x,y) and rearrange to put in
		 * the form of the quadratic formula: ax^2 + bx + c = 0
		 * double m = sin(track->getPhi()) / cos(track->getPhi());
		 */
		double m = sin(angle) / cos(angle);
		q = y0 - m * x0;
		a = _A + _B * _B * m * m;
		b = 2 * _B * m * q + _C + _D * m;
		c = _B * q * q + _D * q + _E;

		discr = b*b - 4*a*c;

		/* There are no intersections */
		if (discr < 0)
			return 0;

		/* There is one intersection (ie on the surface) */
		else if (discr == 0) {
			xcurr = -b / (2*a);
			ycurr = y0 + m * (points[0].getX() - x0);
			points[num].setCoords(xcurr, ycurr);
			if (angle < M_PI && ycurr > y0)
				num++;
			else if (angle > M_PI && ycurr < y0)
				num++;
			return num;
		}

		/* There are two intersections */
		else {
			xcurr = (-b + sqrt(discr)) / (2*a);
			ycurr = y0 + m * (xcurr - x0);
			points[num].setCoords(xcurr, ycurr);
			if (angle < M_PI && ycurr > y0) {
				num++;
			}
			else if (angle > M_PI && ycurr < y0) {
				num++;
			}

			xcurr = (-b - sqrt(discr)) / (2*a);
			ycurr = y0 + m * (xcurr - x0);
			points[num].setCoords(xcurr, ycurr);
			if (angle < M_PI && ycurr > y0) {
				num++;
			}
			else if (angle > M_PI && ycurr < y0) {
				num++;
			}

			return num;
		}
	}
}


/**
 * Converts this Plane's attributes to a character array
 * @param a character array of this plane's attributes
 */
std::string Circle::toString() {
	std::stringstream string;

	string << "Surface id = " << _id << ", type = CIRCLE " << ", A = "
			<< _A << ", B = " << _B << ", C = " << _C << ", D = " << _D
			<< ", E = " << _E;

	return string.str();
}

double Circle::getXMin(){
	log_printf(ERROR, "Circle::getXMin not implemented");
	return std::numeric_limits<double>::infinity();
}

double Circle::getXMax(){
	log_printf(ERROR, "Circle::getXMax not implemented");
	return std::numeric_limits<double>::infinity();
}

double Circle::getYMin(){
	log_printf(ERROR, "Circle::getYMin not implemented");
	return std::numeric_limits<double>::infinity();
}

double Circle::getYMax(){
	log_printf(ERROR, "Circle::getYMax not implemented");
	return std::numeric_limits<double>::infinity();
}


