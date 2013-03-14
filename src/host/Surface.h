/*
 * Surface.h
 *
 *  Created on: Jan 9, 2012
 *      Author: wbinventor
 */


#ifndef SURFACE_H_
#define SURFACE_H_

#include <limits>
#include "Cell.h"
#include "LocalCoords.h"


/* Define for compiler */
class LocalCoords;
class Cell;

/**
 * Surface types
 */
enum surfaceType {
	PLANE,
	CIRCLE,
	XPLANE,
	YPLANE,
	ZPLANE,
	QUADRATIC
};


/**
 * Surface boundary types
 */
enum boundaryType {
	BOUNDARY_NONE,
	REFLECTIVE,
	VACUUM
};



/**
 * Represents a 2-dimensional quadratics surface
 */
class Surface {

protected:

	static short int _n;		  /* Counts the number of surfaces */
	short int _uid;		  /* monotonically increasing id based on n */
	short int _id;
	surfaceType _type;
	boundaryType _boundary;

public:

	Surface(const short int id, const surfaceType type,
		const boundaryType boundary);
	virtual ~Surface();
	short int getUid() const;
	short int getId() const;
	surfaceType getType() const;
	boundaryType getBoundary();
	virtual double evaluate(const Point* point) const =0;
	virtual int intersection(Point* point, double angle, Point* points) =0;
	virtual std::string toString() =0;
	virtual double getXMin() =0;
	virtual double getXMax() =0;
	virtual double getYMin() =0;
	virtual double getYMax() =0;
	bool onSurface(Point* point);
	bool onSurface(LocalCoords* coord);
	double getMinDistance(Point* point, double angle, Point* intersection);

};


/**
 * Represents a plane in 2D as a Surface subclass
 */
class Plane: public Surface {

protected:

	double _A, _B, _C, _D;
	friend class Surface;
	friend class Circle;

public:

	Plane(const short int id, const boundaryType boundary, const double A,
							const double B, const double C, const double D);
	double evaluate(const Point* point) const;
	int intersection(Point* point, double angle, Point* points);
	std::string toString();
	virtual double getXMin();
	virtual double getXMax();
	virtual double getYMin();
	virtual double getYMax();

};


/**
 * Represents a plane parallel to the x-axis as a Plane subclass
 */
class XPlane: public Plane {

private:

public:

	XPlane(const short int id, const boundaryType boundary, const double D);
	std::string toString();
	virtual double getXMin();
	virtual double getXMax();
	virtual double getYMin();
	virtual double getYMax();

};


/**
 * Represents a plane parallel to the y-axis as a Plane subclass
 */
class YPlane: public Plane {

private:

public:

	YPlane(const short int id, const boundaryType boundary, const double D);
	std::string toString();
	virtual double getXMin();
	virtual double getXMax();
	virtual double getYMin();
	virtual double getYMax();

};


/**
 * Represents a plane parallel to the z-axis as a Plane subclass
 */
class ZPlane: public Plane {

private:

public:

	ZPlane(const short int id, const boundaryType boundary, const double D);
	std::string toString();
	virtual double getXMin();
	virtual double getXMax();
	virtual double getYMin();
	virtual double getYMax();

};


/**
 * Represents a circle as a Surface subclass
 */
class Circle: public Surface {

private:

	Point center;
	double _radius;
	double _A, _B, _C, _D, _E;
	friend class Surface;
	friend class Plane;

public:

	Circle(const short int id, const boundaryType boundary, const double x,
				const double y, const double radius);
	double evaluate(const Point* point) const;
	int intersection(Point* point, double angle, Point* points);
	std::string toString();
	virtual double getXMin();
	virtual double getXMax();
	virtual double getYMin();
	virtual double getYMax();
	double getRadius();

};


/**
 * Finds the minimum distance from a point with a given trajectory defined
 * by an angle to this surface. If the trajectory will not intersect the
 * surface, returns INFINITY
 * @param point a pointer to the point of interest
 * @param angle the angle defining the trajectory in radians
 * @param intersection a pointer to a point for storing the intersection
 * @return the minimum distance
 */
inline double Surface::getMinDistance(Point* point, double angle,
									Point* intersection) {

	/* Point array for intersections with this surface */
	Point intersections[2];

	/* Find the intersection point(s) */
	int num_inters = this->intersection(point, angle, intersections);
	double distance = INFINITY;

	/* If there is one intersection point */
	if (num_inters == 1) {
		distance = intersections[0].distance(point);
		intersection->setX(intersections[0].getX());
		intersection->setY(intersections[0].getY());
	}

	/* If there are two intersection points */
	else if (num_inters == 2) {
		double dist1 = intersections[0].distance(point);
		double dist2 = intersections[1].distance(point);

		/* Determine which intersection point is nearest */
		if (dist1 < dist2) {
			distance = dist1;	double getRadius();

			intersection->setX(intersections[0].getX());
			intersection->setY(intersections[0].getY());
		}
		else {
			distance = dist2;
			intersection->setX(intersections[1].getX());
			intersection->setY(intersections[1].getY());
		}
	}

	return distance;
}


/**
 * Evaluate a point using the plane's quadratic surface equation
 * @param point a pointer to the point of interest
 * @return the value of point in the equation
 */
inline double Plane::evaluate(const Point* point) const {
	double x = point->getX();
	double y = point->getY();
	//TODO: does not support z-planes
	return (_A * x + _B * y + _D);
}


/**
 * Return the radius of the circle
 * @return the radius of the circle
 */
inline double Circle::getRadius() {
	return this->_radius;
}


/**
 * Evaluate a point using the circle's quadratic surface equation
 * @param point a pointer to the point of interest
 * @return the value of point in the equation
 */
inline double Circle::evaluate(const Point* point) const {
	double x = point->getX();
	double y = point->getY();
	return (_A * x * x + _B * y * y + _C * x + _D * y + _E);
}


#endif /* SURFACE_H_ */
