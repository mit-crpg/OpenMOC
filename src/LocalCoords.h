/**
 * @file LocalCoords.h
 * @brief The LocalCoords class.
 * @date January 25, 2012
 * @author William Boyd, MIT, Course 22 (wboyd@mit.edu)
 */

#ifndef LOCALCOORDS_H_
#define LOCALCOORDS_H_

#ifdef __cplusplus
#include "Point.h"
#include "Universe.h"
#endif

/**
 * @enum coordType
 * @brief The type of universe level on which the localcoords reside
 */
enum coordType {
    /** A universe level coordinate type */
    UNIV,
    /** A lattice level coordinate type */
    LAT
};


/**
 * @class LocalCoords LocalCoords.h "openmoc/src/host/LocalCoords.h"
 * @brief The localcoords represents a set of local coordinates on some
 *        level of nested universes making up the geometry.
 */
class LocalCoords {

private:
    /** The local coordinate type (UNIV or LAT) */
    coordType _type;
    /** The id of the universe within which this localcoords resides */
    int _universe;
    /** The id of the cell within which this localcoords resides */
    int _cell;
    /** The id of the lattice within which this localcoords resides */
    int _lattice;
    /** The first index of the lattice cell within which this localcoords 
     *  resides */
    int _lattice_x;
    /** The second index of the lattice cell within which this localcoords 
     *  resides */
    int _lattice_y;
    /** A point representing the 2D coordinates of this localcoords */
    Point _coords;
    /** A pointer to the localcoords at the next lower nested universe level */
    LocalCoords* _next;
    /** A pointer to the localcoords at the next higher nested universe level */
    LocalCoords* _prev;

public:
    LocalCoords(double x, double y);
    virtual ~LocalCoords();
    coordType getType();
    int getUniverse() const;
    int getCell() const;
    int getLattice() const;
    int getLatticeX() const;
    int getLatticeY() const;
    double getX() const;
    double getY() const;
    Point* getPoint();
    LocalCoords* getNext() const;
    LocalCoords* getPrev() const;

    void setType(coordType type);
    void setUniverse(int universe);
    void setCell(int cell);
    void setLattice(int lattice);
    void setLatticeX(int lattice_x);
    void setLatticeY(int lattice_y);
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
