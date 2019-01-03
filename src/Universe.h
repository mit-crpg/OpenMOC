/**
 * @file Universe.h
 * @brief The Universe class.
 * @date January 9, 2012
 * @author William Boyd, MIT, Course 22 (wboyd@mit.edu)
 */

#ifndef UNIVERSE_H_
#define UNIVERSE_H_

#ifdef __cplusplus
#ifdef SWIG
#include "Python.h"
#endif
#include "constants.h"
#include "LocalCoords.h"
#include "boundary_type.h"
#include <limits>
#include <map>
#include <vector>
#endif


/* Forward declarations to resolve circular dependencies */
class LocalCoords;
class Cell;
class Surface;
class Material;


int universe_id();
void reset_universe_id();
void maximize_universe_id(int universe_id);


/**
 * @enum universeType
 * @brief The type of universe.
 */
enum universeType{

  /** A simple non-repeating Universe */
  SIMPLE,

  /** A collection of Universes in a rectangular Lattice */
  LATTICE
};


/**
 * @class Universe Universe.h "src/Universe.h"
 * @brief A Universe represents an unbounded space in 3D.
 * @details A Universe contains cell which are bounded subspaces in 3D
 *          which together form the Universe. Universes allow
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

  /** A user-defined name for the Surface */
  char* _name;

  /** The type of Universe (ie, SIMPLE or LATTICE) */
  universeType _type;

  /** A collection of Cell IDs and Cell pointers in this Universe */
  std::map<int, Cell*> _cells;

  /** A boolean representing whether or not this Universe contains a Material
   *  with a non-zero fission cross-section and is fissionable */
  bool _fissionable;

  /** The extrema of the Universe */
  double _min_x;
  double _max_x;
  double _min_y;
  double _max_y;
  double _min_z;
  double _max_z;

  /** A flag for determining if boundaries are up to date */
  bool _boundaries_inspected;

  /** The boundaryTypes of the universe */
  boundaryType _min_x_bound;
  boundaryType _max_x_bound;
  boundaryType _min_y_bound;
  boundaryType _max_y_bound;
  boundaryType _min_z_bound;
  boundaryType _max_z_bound;

public:

  Universe(const int id=-1, const char* name="");
  virtual ~Universe();
  int getUid() const;
  int getId() const;
  char* getName() const;
  universeType getType();
  int getNumCells() const;
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

  Cell* getCell(int cell_id);
  std::map<int, Cell*> getCells() const;
  std::map<int, Cell*> getAllCells();
  std::map<int, Material*> getAllMaterials();
  std::map<int, Universe*> getAllUniverses();
  bool isFissionable();

  void resetBoundaries();
  void calculateBoundaries();
  void setName(const char* name);
  void setType(universeType type);
  void addCell(Cell* cell);
  void removeCell(Cell* cell);

  bool containsPoint(Point* point);
  Cell* findCell(LocalCoords* coords);
  void setFissionability(bool fissionable);
  void subdivideCells(double max_radius=INFINITY);
  void buildNeighbors();

  virtual std::string toString();
  void printString();

  Universe* clone();
};


/**
 * @class Lattice Universe.h "src/Universe.h"
 * @brief Represents a repeating 3D Lattice of Universes.
 */
class Lattice: public Universe {

private:

  /** The number of Lattice cells along the x-axis */
  int _num_x;

  /** The number of Lattice cells along the y-axis */
  int _num_y;

  /** The number of Lattice cells along the z-axis */
  int _num_z;

  /** True if the lattice is non-uniform */
  bool _non_uniform;

  /** The width of each Lattice cell (cm) along the x-axis
      (uniform lattices only) */
  double _width_x;

  /** x-direction dimensions of non-uniform lattice meshes */
  std::vector<double> _widths_x;
  std::vector<double> _accumulate_x;

  /** The width of each Lattice cell (cm) along the y-axis 
      (uniform lattices only) */
  double _width_y;

  /** y-direction dimensions of non-uniform lattice meshes */
  std::vector<double> _widths_y;
  std::vector<double> _accumulate_y;

  /** The width of each Lattice cell (cm) along the z-axis 
      (uniform lattices only) */
  double _width_z;

  /** z-direction dimensions of non-uniform lattice meshes */
  std::vector<double> _widths_z;
  std::vector<double> _accumulate_z;

  /** The coordinates of the offset for the Universe */
  Point _offset;

  /** A container of Universes */
  std::vector< std::vector< std::vector< std::pair<int, Universe*> > > >
      _universes;

public:

  Lattice(const int id=-1, const char* name="");
  virtual ~Lattice();

  void setOffset(double x, double y, double z=0.0);
  Point* getOffset();
  int getNumX() const;
  int getNumY() const;
  int getNumZ() const;
  double getWidthX() const;
  double getWidthY() const;
  double getWidthZ() const;
  bool getNonUniform() const;
  const std::vector<double>& getWidthsX() const;
  const std::vector<double>& getWidthsY() const;
  const std::vector<double>& getWidthsZ() const;
  const std::vector<double>& getAccumulateX() const;
  const std::vector<double>& getAccumulateY() const;
  const std::vector<double>& getAccumulateZ() const;
  double getMinX();
  double getMaxX();
  double getMinY();
  double getMaxY();
  double getMinZ();
  double getMaxZ();

  Universe* getUniverse(int lat_x, int lat_y, int lat_z=0) const;
  std::vector< std::vector< std::vector< std::pair<int, Universe*> > > >*
      getUniverses();
  std::map<int, Universe*> getUniqueUniverses();
  std::map<int, double> getUniqueRadius(std::map<int, Universe*> unique_universes);
  std::map<int, Cell*> getAllCells();
  std::map<int, Universe*> getAllUniverses();

  void setNumX(int num_x);
  void setNumY(int num_y);
  void setNumZ(int num_z);
  void setWidth(double width_x, double width_y,
                double width_z=std::numeric_limits<double>::infinity());
  void setNonUniform(bool non_uniform);
  void setWidthsX(std::vector<double> widthsx);
  void setWidthsY(std::vector<double> widthsy);
  void setWidthsZ(std::vector<double> widthsz);
  void setAccumulateX(std::vector<double> accumulatex);
  void setAccumulateY(std::vector<double> accumulatey);
  void setAccumulateZ(std::vector<double> accumulatez);
  void setUniverses(int num_z, int num_y, int num_x, Universe** universes);
  void updateUniverse(int lat_x, int lat_y, int lat_z, Universe* universe);
  void removeUniverse(Universe* universe);
  void subdivideCells(double max_radius=INFINITY);
  void buildNeighbors();

  bool containsPoint(Point* point);
  Cell* findCell(LocalCoords* coords);
  double minSurfaceDist(Point* point, double azim, double polar=M_PI/2.0);

  int getLatX(Point* point);
  int getLatY(Point* point);
  int getLatZ(Point* point);

  int getLatticeCell(Point* point);
  int getLatticeSurface(int cell, Point* point);
  int getLatticeSurfaceOTF(int cell, double z, int surface_2D);

  std::string toString();
  void printString();

  /* Set XYZ widths of non-uniform meshes */
  void setWidths(std::vector<double> widths_x, std::vector<double> widths_y, 
                 std::vector<double> widths_z);
  void computeSizes();

  /* For debug use */
  void printLatticeSizes();

};

/**
 * @brief A helper struct for the Universe::findCell() method.
 * @details This is used to insert a Universe's Cells to the back of a vector
 *          of neighbor Cells in Universe::findCell() routine. This works in
 *          symbiosis with the pair_second method template defined below.
 */
template<typename tPair>
struct second_t {
  typename tPair::second_type operator()(const tPair& p) const {
    return p.second;
  }
};


/**
 * @brief A helper routine for the Universe::findCell() method.
 * @details This is used to insert a Universe's Cells to the back of a vector
 *          of neighbor Cells in Universe::findCell() routine. This works in
 *          symbiosis with the second_t struct template defined above.
 * @param map a std::map iterator
 * @return the second element in the iterator (e.g., map value)
 */
template<typename tMap>
second_t<typename tMap::value_type> pair_second(const tMap& map) {
  return second_t<typename tMap::value_type>();
}

#endif /* UNIVERSE_H_ */
