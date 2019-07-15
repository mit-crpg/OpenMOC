/**
 * @file LocalCoords.h
 * @brief The LocalCoords class.
 * @date January 25, 2012
 * @author William Boyd, MIT, Course 22 (wboyd@mit.edu)
 */

#ifndef LOCALCOORDS_H_
#define LOCALCOORDS_H_

#ifdef __cplusplus
#ifdef SWIG
#include "Python.h"
#endif
#include "Point.h"
#include "Universe.h"
#include "Cell.h"
#include "constants.h"
#endif

/* Forward declarations to resolve circular dependencies */
class Universe;
class Lattice;
class Cell;

/**
 * @enum coordType
 * @brief The type of Universe level on which the LocalCoords reside
 */
enum coordType {
  /** A Universe level coordinate type */
  UNIV,

  /** A Lattice level coordinate type */
  LAT
};


/**
 * @class LocalCoords LocalCoords.h "openmoc/src/host/LocalCoords.h"
 * @brief The LocalCoords represents a set of local coordinates on some
 *        level of nested Universes making up the geometry.
 *        _next and _prev allow for use of LocalCoords as a linked list
 *        but _next_array can also be used to access coordinates.
 */
class LocalCoords {

private:
  /** The local coordinate type (UNIV or LAT) */
  coordType _type;

  /** The Universe within which this LocalCoords resides */
  Universe* _universe;

  /** The Cell within which this LocalCoords resides */
  Cell* _cell;

  /** The Lattice within which this LocalCoords resides */
  Lattice* _lattice;

  /** The first index of the Lattice cell within which this LocalCoords
   *  resides */
  int _lattice_x;

  /** The second index of the Lattice cell within which this LocalCoords
   *  resides */
  int _lattice_y;

  /** The third index of the Lattice cell within which this LocalCoords
   *  resides */
  int _lattice_z;

  /** A Point representing the 3D coordinates of this LocalCoords */
  Point _coords;

  /** The direction angle in radians with respect to the x-axis */
  double _phi;

  /** The direction angle in radians with respect to the z-axis */
  double _polar;

  /** A pointer to the LocalCoords at the next lower nested Universe level */
  LocalCoords* _next;

  /** A pointer to the LocalCoords at the next higher nested Universe level */
  LocalCoords* _prev;

  /** An array that contains pointers to all the next LocalCoords */
  LocalCoords* _next_array;

  /** Position in the _next_array of coordinates */
  int _position;

  /** Size of the _next_array of coordinates */
  int _array_size;

  /** An integer to differentiate otherwise matching coordinate FSR keys */
  int _version_num;

  void setArrayPosition(LocalCoords* array, int position, int array_size);

public:
  LocalCoords(double x=0.0, double y=0.0, double z=0.0, bool first=false);
  virtual ~LocalCoords();
  coordType getType();
  Universe* getUniverse() const;
  Cell* getCell() const;
  Lattice* getLattice() const;
  int getLatticeX() const;
  int getLatticeY() const;
  int getLatticeZ() const;
  double getX() const;
  double getY() const;
  double getZ() const;
  double getPhi() const;
  double getPolar() const;
  Point* getPoint();
  LocalCoords* getNext() const;
  LocalCoords* getNextCreate(double x, double y, double z);
  LocalCoords* getPrev() const;
  int getVersionNum();
  int getPosition();

  void setType(coordType type);
  void setUniverse(Universe* universe);
  void setCell(Cell* cell);
  void setLattice(Lattice* lattice);
  void setLatticeX(int lattice_x);
  void setLatticeY(int lattice_y);
  void setLatticeZ(int lattice_z);
  void setX(double x);
  void setY(double y);
  void setZ(double z);
  void setPhi(double phi);
  void setPolar(double polar);
  void setNext(LocalCoords *next);
  void setPrev(LocalCoords* coords);
  void setVersionNum(int version_num);

  LocalCoords* getLowestLevel();
  LocalCoords* getHighestLevel();
  void adjustCoords(double delta_x, double delta_y, double delta_z=0.0);
  void updateMostLocal(Point* point);
  void prune();
  void deleteArray();
  void copyCoords(LocalCoords* coords);
  std::string toString();
  void detectLoop();
};


#endif /* LOCALCOORDS_H_ */
