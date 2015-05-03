/**
 * @file Cell.h
 * @brief The Cell class.
 * @date January 18, 2012
 * @author William Boyd, MIT, Course 22 (wboyd@mit.edu)
 */


#ifndef CELL_H_
#define CELL_H_

#ifdef __cplusplus
#include "Python.h"
#include "Material.h"
#include "Surface.h"
#include "Point.h"
#include <limits>
#include <map>
#include <vector>
#endif

/* Forward declarations to resolve circular dependencies */
class Universe;
class Surface;

int cell_id();
void reset_cell_id();


/**
 * @struct surface_halfspace
 * @brief A surface_halfspace represents a surface pointer with associated
 *        halfspace.
 */
struct surface_halfspace {

  /** A pointer to the Surface object */
  Surface* _surface;

  /** The halfspace associated with this surface */
  int _halfspace;

};



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

  /** The ID for the Universe within which this cell resides */
  int _universe;

  /** Map of bounding Surface IDs with pointers and halfspaces (+/-1) */
  std::map<int, surface_halfspace> _surfaces;

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

  /** The boundary condition at the minimum reachable x-coordinate */
  boundaryType _min_x_bc;

  /** The boundary condition at the maximum reachable x-coordinate */
  boundaryType _max_x_bc;

  /** The boundary condition at the minimum reachable y-coordinate */
  boundaryType _min_y_bc;

  /** The boundary condition at the maximum reachable y-coordinate */
  boundaryType _max_y_bc;

  /** The boundary condition at the minimum reachable z-coordinate */
  boundaryType _min_z_bc;

  /** The boundary condition at the maximum reachable z-coordinate */
  boundaryType _max_z_bc;

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
  boundaryType getMinXBoundaryType();
  boundaryType getMaxXBoundaryType();
  boundaryType getMinYBoundaryType();
  boundaryType getMaxYBoundaryType();
  boundaryType getMinZBoundaryType();
  boundaryType getMaxZBoundaryType();
  int getNumSurfaces() const;
  std::map<int, surface_halfspace> getSurfaces() const;

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
  std::map<int, Cell*> getAllCells();
  std::map<int, Universe*> getAllUniverses();

  void setMaterial(Material* material);
  void setNumRings(int num_rings);
  void setNumSectors(int num_sectors);

  CellBasic* clone();
  std::vector<CellBasic*> subdivideCell();

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
  std::map<int, Cell*> getAllCells();
  std::map<int, Universe*> getAllUniverses();

  void setFill(Universe* universe_fill);

  std::string toString();
  void printString();
};

#endif /* CELL_H_ */
