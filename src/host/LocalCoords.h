/*
 * LocalCoords.h
 *
 *  Created on: Jan 25, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#ifndef LOCALCOORDS_H_
#define LOCALCOORDS_H_

#include "Point.h"
#include "Universe.h"


/* Type represents whether a localcoords is in a simple
 * universe or lattice */
enum coordType {
	UNIV,
	LAT
};


class LocalCoords {

private:

	coordType _type;
	short int _universe;
	short int _cell;
	short int _lattice;
	short int _lattice_x;
	short int _lattice_y;
	Point _coords;
	LocalCoords* _next;
	LocalCoords* _prev;

public:

	LocalCoords(double x, double y);
	virtual ~LocalCoords();
	coordType getType();
    short int getUniverse() const;
    short int getCell() const;
    short int getLattice() const;
    short int getLatticeX() const;
    short int getLatticeY() const;
    double getX() const;
    double getY() const;
    Point* getPoint();
    LocalCoords* getNext() const;
    LocalCoords* getPrev() const;

    void setType(coordType type);
    void setUniverse(short int universe);
    void setCell(short int cell);
    void setLattice(short int lattice);
    void setLatticeX(short int lattice_x);
    void setLatticeY(short int lattice_y);
    void setX(double x);
    void setY(double y);
    void setNext(LocalCoords *next);
    void setPrev(LocalCoords* coords);

    LocalCoords* getLowestLevel();
    void adjustCoords(double delta_x, double delta_y);
    void updateMostLocal(Point* point);
    void prune();
    void copyCoords(LocalCoords* coords);
    std::string toString();

};

#endif /* LOCALCOORDS_H_ */
