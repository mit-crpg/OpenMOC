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
#include "Region.h"
#include "Surface.h"
#include "Point.h"
#include <limits>
#include <string>
#endif

/* Forward declarations to resolve circular dependencies */
class Universe;
class Surface;
class Region;
class Halfspace;

int cell_id();
void reset_cell_id();
void maximize_cell_id(int cell_id);


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

  /** A pointer to the Region bounding the Cell */
  Region* _region;

  /** A pointer to the Region currently examined while transfering geometry */
  Region* _current_region;

  /** The volume / area of the Cell computed from overlapping segments */
  double _volume;

  /** The total number of instances of this Cell in the Geometry */
  int _num_instances;

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

  /** A parent Cell if cloned by another Cell */
  Cell* _parent;

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
  Region* getRegion();
  double getVolume();
  int getNumInstances();
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
  boundaryType getMinZBoundaryType();
  boundaryType getMaxZBoundaryType();
  int getNumSurfaces() const;
  std::map<int, Halfspace*> getSurfaces() const;
  std::vector<Cell*> getNeighbors() const;
  bool hasParent();
  Cell* getParent();
  Cell* getOldestAncestor();

  int getNumZCylinders();

  std::map<int, Cell*> getAllCells();
  std::map<int, Universe*> getAllUniverses();

  void setName(const char* name);
  void setFill(Material* fill);
  void setFill(Universe* fill);
  void setRegion(Region* region);
  void setVolume(double volume);
  void incrementVolume(double volume);
  void setNumInstances(int num_instances);
  void incrementNumInstances();
  void setRotation(double* rotation, int num_axes, std::string units="degrees");
  void setTranslation(double* translation, int num_axes);
  void setNumRings(int num_rings, double inner_radius=-1);
  void setNumSectors(int num_sectors);
  void setParent(Cell* parent);
  void addSurface(int halfspace, Surface* surface);
  void addSurfaceInRegion(int halfspace, Surface* surface);
  void addLogicalNode(int region_type);
  void goUpOneRegionLogical();
  void removeSurface(Surface* surface);
  void addNeighborCell(Cell* cell);

  bool isFissionable();
  bool containsPoint(Point* point);
  bool containsCoords(LocalCoords* coords);
  double minSurfaceDist(Point* point, double azim, double polar);
  double minSurfaceDist(LocalCoords* coords);

  Cell* clone(bool clone_region=true);
  void subdivideCell(double max_radius);
  void buildNeighbors();

  std::string toString();
  void printString();
};


#endif /* CELL_H_ */
