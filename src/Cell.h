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

  /** A cell filled by a Material */
  MATERIAL,

  /** A cell filled by a Universe */
  FILL
};


/**
 * @class Cell Cell.h "src/Cell.h"
 * @brief Represents a Cell inside of a Universe.
 */
class Cell {

private:

  /** A static counter for the number of Cells */
  static int _n;

  /** A monotonically increasing unique ID for each Cell created */
  int _uid;

  /** A user-defined ID for each Cell created */
  int _id;

  /** A user-defined name for the Surface */
  char* _name;

  /** The type of Cell (ie MATERIAL or FILL) */
  cellType _cell_type;

  /** A pointer to the Material or Universe filling this Cell */
  void* _fill;

  /** The number of rings sub-dividing this Cell */
  int _num_rings;

  /** The number of sectors sub-dividing this Cell */
  int _num_sectors;

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

  /** A container of all Cell clones created for rings */
  std::vector<Cell*> _rings;

  /** A container of all Cell clones created for angular sectors */
  std::vector<Cell*> _sectors;

  /** A container of all Cell clones created for rings and sectors */
  std::vector<Cell*> _subcells;

  void ringify();
  void sectorize();

public:
  Cell(int id=0, const char* name="");
  virtual ~Cell();
  int getUid() const;
  int getId() const;
  char* getName() const;
  cellType getType() const;
  Material* getFillMaterial();
  Universe* getFillUniverse();
  int getNumRings();
  int getNumSectors();
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

  std::map<int, Cell*> getAllCells();
  std::map<int, Universe*> getAllUniverses();

  void setName(const char* name);
  void setFill(Material* fill);
  void setFill(Universe* fill);
  void setNumRings(int num_rings);
  void setNumSectors(int num_sectors);
  void addSurface(int halfspace, Surface* surface);
  void removeSurface(Surface* surface);
  void findBoundingBox();
  bool cellContainsPoint(Point* point);
  bool cellContainsCoords(LocalCoords* coords);
  double minSurfaceDist(Point* point, double angle, Point* min_intersection);

  Cell* clone();
  std::vector<Cell*> subdivideCell();

  std::string toString();
  void printString();
};


#endif /* CELL_H_ */
