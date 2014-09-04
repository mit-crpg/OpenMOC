/**
 * @file Universe.h
 * @brief The Universe class.
 * @date January 9, 2012
 * @author William Boyd, MIT, Course 22 (wboyd@mit.edu)
 */

#ifndef UNIVERSE_H_
#define UNIVERSE_H_

#ifdef __cplusplus
#include <map>
#include <vector>
#include "Cell.h"
#include "LocalCoords.h"
#endif


/** Error threshold for determining how close to the boundary of a Lattice cell
 * a Point needs to be to be considered on it */
#define ON_LATTICE_CELL_THRESH 1E-12


/** Distance a Point is moved to cross over a Surface into a new Cell during
 * Track segmentation */
#define TINY_MOVE 1E-10


class LocalCoords;
class Cell;
class CellFill;
class CellBasic;


int universe_id();


/**
 * @enum universeType
 * @brief The type of universe
 */
enum universeType{

  /** A simple non-repeating Universe */
  SIMPLE,

  /** A collection of Universes in a rectangular Lattice */
  LATTICE
};


/**
 * @class Universe Universe.h "src/Universe.h"
 * @brief A Universe represents an unbounded space in the 2D xy-plane.
 * @details A Universe contains cell which are bounded subspaces in the 2D
 *          xy-plane and which together form the Universe. Universes allow
 *          for complex, repeating (i.e. lattices) geometries to be simply
 *          represented with as few data structures as possible.
 */
class Universe {

protected:

  /** A static counter for the number of Universes */
  static int _n;

  /** A monotonically increasing unique ID for each Universe created */
  int _uid;

  /** A user-defined id for each Universe created */
  int _id;

  /** The type of Universe (ie, SIMPLE or LATTICE) */
  universeType _type;

  /** A collection of Cell IDs and Cell pointers */
  std::map<int, Cell*> _cells;

  /** A boolean representing whether or not this Universe contains a Material
   *  with a non-zero fission cross-section and is fissionable */
  bool _fissionable;

public:

  Universe(const int id);
  virtual ~Universe();

  void addCell(Cell* cell);

  Cell* getCell(int cell_id);
  CellFill* getCellFill(int cell_id);
  CellBasic* getCellBasic(int cell_id);
  std::map<int, Cell*> getCells() const;
  int getUid() const;
  int getId() const;
  universeType getType();
  int getNumCells() const;
  std::vector<int> getMaterialIds();
  std::vector<int> getNestedUniverseIds();
  void getCellIds(int* cell_ids, int num_cells);
  bool isFissionable();

  void setType(universeType type);

  void setFissionability(bool fissionable);

  Cell* findCell(LocalCoords* coords, std::map<int, Universe*> universes);
  double minSurfaceDist(Point* point, double angle);
  void subdivideCells();
  std::string toString();
  void printString();

  Universe* clone();
};


/**
 * @class Lattice Universe.h "src/Universe.h"
 * @brief Represents a repeating 2D Lattice of Universes.
 */
class Lattice: public Universe {

private:

  /** The number of Lattice cells along the x-axis */
  int _num_x;

  /** The number of Lattice cells along the y-axis */
  int _num_y;

  /** The width of each Lattice cell (cm) along the x-axis */
  double _width_x;

  /** The width of each Lattice cell (cm) along the y-axis */
  double _width_y;

  /** The coordinates of the offset for the Universe */
  Point _offset;

  /** A container of Universes ? */
  std::vector< std::vector< std::pair<int, Universe*> > > _universes;

public:

  Lattice(const int id, const double width_x, const double width_y);
  virtual ~Lattice();

  void setOffset(double x, double y);
  Point* getOffset();
  int getNumX() const;
  int getNumY() const;
  void setNumX(int num_x);
  void setNumY(int num_y);
  std::vector< std::vector< std::pair<int, Universe*> > >
                                           getUniverses() const;
  Universe* getUniverse(int lattice_x, int lattice_y) const;
  double getWidthX() const;
  double getWidthY() const;
  std::vector<int> getNestedUniverseIds();

  int getLatX(Point* point);
  int getLatY(Point* point);

  int getLatticeCell(Point* point);
  int getLatticeSurface(int cell, Point* point);

  void setLatticeCells(int num_x, int num_y, int* universes);
  void setUniversePointer(Universe* universe);

  bool withinBounds(Point* point);
  Cell* findCell(LocalCoords* coords, std::map<int, Universe*> universes);
  double minSurfaceDist(Point* point, double angle);
  std::string toString();
  void printString();
};

#endif /* UNIVERSE_H_ */

