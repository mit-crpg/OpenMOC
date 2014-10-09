/**
 * @file Cell.h
 * @brief The Cell class.
 * @date January 18, 2012
 * @author William Boyd, MIT, Course 22 (wboyd@mit.edu)
 */


#ifndef CELL_H_
#define CELL_H_

#ifdef __cplusplus
#include <limits>
#include <map>
#include <vector>
#include "Material.h"
#include "Surface.h"
#include "Point.h"
#endif

class Universe;
class Surface;
class LocalCoords;

int cell_id();
void reset_cell_id();


/**
 * @enum cellType
 * @brief The type of cell.
*/
enum cellType {

  /** A cell filled by a material */
  MATERIAL,

  /** A cell filled by a universe */
  FILL
};


/**
 * @class Cell Cell.h "src/Cell.h"
 * @brief Represents a Cell inside of a Universe.
 */
class Cell {

protected:

  /** A static counter for the number of Cells */
  static int _n;

  /** A static counter for the number of times this Cell has been cloned */
  static int _num_clones;

  /** A monotonically increasing unique ID for each Cell created */
  int _uid;

  /** A user-defined ID for each Cell created */
  int _id;

  /** A user-defined name for the Surface */
  char* _name;

  /** The type of Cell (ie MATERIAL or FILL) */
  cellType _cell_type;

  /** Map of bounding Surface pointers to halfspaces (+/-1) */
  std::map<Surface*, int> _surfaces;

  /** A boolean representing whether or not this Cell contains a Material
   *  with a non-zero fission cross-section and is fissionable */
  bool _fissionable;

  /** The minimum reachable x-coordinate within the Cell */
  double _min_x;

  /** The maximum reachable x-coordinate within the Cell */
  double _max_x;

  /** The minimum reachable y-coordinate within the Cell */
  double _min_y;

  /** The maximum reachable y-coordinate within the Cell */
  double _max_y;

  /** The minimum reachable z-coordinate within the Cell */
  double _min_z;

  /** The maximum reachable z-coordinate within the Cell */
  double _max_z;

public:
  Cell();
  Cell(int id=0, const char* name="");
  virtual ~Cell();
  int getUid() const;
  int getId() const;
  char* getName() const;
  cellType getType() const;
  double getMinX();
  double getMaxX();
  double getMinY();
  double getMaxY();
  double getMinZ();
  double getMaxZ();
  int getNumSurfaces() const;
  std::map<Surface*, int> getSurfaces() const;
  bool isFissionable();

  /**
   * @brief Return the number of flat source regions in this Cell.
   * @details This method is used when the Geometry recursively constructs
   *          flat source regions.
   * @return the number of FSRs in this Cell
   */
  virtual int getNumFSRs() =0;

  /**
   * @brief Returns the std::map of Cell IDs and Cell pointers within any
   *        nested Universes filling this Cell.
   * @return std::map of Cell IDs and pointers
   */
  virtual std::map<int, Cell*> getAllCells() =0;

  /**
   * @brief Returns the std::map of Universe IDs and Universe pointers within
   *        any nested Universes filling this Universe.
   * @return std::map of Universe IDs and pointers
   */
  virtual std::map<int, Universe*> getAllUniverses() =0;

  void setName(const char* name);
  void addSurface(int halfspace, Surface* surface);
  void removeSurface(Surface* surface);
  void findBoundingBox();

  /**
   * @brief Computes whether or not this Cell contains a fissionable Material.
   */
  virtual void computeFissionability() =0;

  bool cellContainsPoint(Point* point);
  bool cellContainsCoords(LocalCoords* coords);
  double minSurfaceDist(Point* point, double angle, Point* min_intersection);
  /**
   * @brief Convert this CellFill's attributes to a string format.
   * @return a character array of this Cell's attributes
   */
  virtual std::string toString() =0;

  /**
   * @brief Prints a string representation of all of the Cells's objects to
   *        the console.
   */
  virtual void printString() =0;
};


/**
 * @class CellBasic Cell.h "src/Cell.h"
 * @brief Represents a Cell filled with a Material.
 */
class CellBasic: public Cell {

private:

  /** A pointer to the Material filling this Cell */
  Material* _material;

  /** The number of rings sub-dividing this Cell */
  int _num_rings;

  /** The number of sectors sub-dividing this Cell */
  int _num_sectors;

  /** A container of all CellBasic clones created for rings */
  std::vector<CellBasic*> _rings;

  /** A container of all CellBasic clones created for angular sectors */
  std::vector<CellBasic*> _sectors;

  /** A container of all CellBasic clones created for rings and sectors */
  std::vector<CellBasic*> _subcells;

  void ringify();
  void sectorize();

public:
  CellBasic(int id=0, const char* name="", int rings=0, int sectors=0);

  Material* getMaterial();
  int getNumRings();
  int getNumSectors();
  int getNumFSRs();
  std::map<int, Cell*> getAllCells();
  std::map<int, Universe*> getAllUniverses();

  void setMaterial(Material* material);
  void setNumRings(int num_rings);
  void setNumSectors(int num_sectors);

  CellBasic* clone();
  std::vector<CellBasic*> subdivideCell();
  void computeFissionability();

  std::string toString();
  void printString();
};


/**
 * @class CellFill Cell.h "src/Cell.h"
 * @brief Represents a Cell filled with a Universe.
 */
class CellFill: public Cell {

private:

  /** The pointer to the Universe filling this Cell */
  Universe* _fill;

public:
  CellFill(int id=0, const char* name="");

  Universe* getFill() const;
  int getNumFSRs();
  std::map<int, Cell*> getAllCells();
  std::map<int, Universe*> getAllUniverses();

  void setFill(Universe* universe_fill);

  void computeFissionability();

  std::string toString();
  void printString();
};

#endif /* CELL_H_ */
