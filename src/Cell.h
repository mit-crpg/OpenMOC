/**
 * @file Cell.h
 * @brief The Cell class.
 * @date January 18, 2012
 * @author William Boyd, MIT, Course 22 (wboyd@mit.edu)
 */


#ifndef CELL_H_
#define CELL_H_

#ifdef __cplusplus
#ifdef SWIG
#include "Python.h"
#endif
#include "Material.h"
#include "Surface.h"
#include "Point.h"
#include <limits>
#include <string>
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
  FILL,

  /** A cell not yet filled by anything */
  UNFILLED
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

  /* A boolean indicating whether to cell is rotated */
  bool _rotated;

  /* An array with angles in degrees for rotations about x, y, and z */
  double _rotation[3];
  
  /* A rotation matrix defined in terms of the rotation angles */
  double _rotation_matrix[9];

  /* A boolean indicating whether to cell is translated */
  bool _translated;

  /* An array with translations in x, y and z */
  double _translation[3];

  /** The number of rings sub-dividing this Cell */
  int _num_rings;

  /** The number of sectors sub-dividing this Cell */
  int _num_sectors;

  /** Map of bounding Surface IDs with pointers and halfspaces (+/-1) */
  std::map<int, surface_halfspace*> _surfaces;

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

  /* Vector of neighboring Cells */
  std::vector<Cell*> _neighbors;

  void ringify(std::vector<Cell*>& subcells, double max_radius);
  void sectorize(std::vector<Cell*>& subcells);

public:
  Cell(int id=0, const char* name="");
  virtual ~Cell();
  int getUid() const;
  int getId() const;
  char* getName() const;
  cellType getType() const;
  Material* getFillMaterial();
  Universe* getFillUniverse();
  bool isRotated();
  bool isTranslated();
  double getPhi(std::string units="degrees");
  double getTheta(std::string units="degrees");
  double getPsi(std::string units="degrees");
  double* getRotationMatrix();
  double* getTranslation();
  void retrieveRotation(double* rotations, int num_axes,
			std::string units="degrees");
  void retrieveTranslation(double* translations, int num_axes);
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
  int getNumSurfaces() const;
  std::map<int, surface_halfspace*> getSurfaces() const;
  std::vector<Cell*> getNeighbors() const;

  std::map<int, Cell*> getAllCells();
  std::map<int, Universe*> getAllUniverses();

  void setName(const char* name);
  void setFill(Material* fill);
  void setFill(Universe* fill);
  void setRotation(double* rotation, int num_axes, std::string units="degrees");
  void setTranslation(double* translation, int num_axes);
  void setNumRings(int num_rings);
  void setNumSectors(int num_sectors);
  void addSurface(int halfspace, Surface* surface);
  void removeSurface(Surface* surface);
  void addNeighborCell(Cell* cell);

  void findBoundingBox();
  bool containsPoint(Point* point);
  bool containsCoords(LocalCoords* coords);
  double minSurfaceDist(LocalCoords* coords);

  Cell* clone();
  void subdivideCell(double max_radius);
  void buildNeighbors();

  std::string toString();
  void printString();
};


#endif /* CELL_H_ */
