/**
 * @file Geometry.h
 * @brief The Geometry class.
 * @date January 19, 2012
 * @author William Boyd, MIT, Course 22 (wboyd@mit.edu)
 */

#ifndef GEOMETRY_H_
#define GEOMETRY_H_

#ifdef __cplusplus
#ifdef SWIG
#include "Python.h"
#endif
#include "Cmfd.h"
#include <limits>
#include <sys/types.h>
#include <sys/stat.h>
#include <sstream>
#include <string>
#include <omp.h>
#include <functional>
#include "ParallelHashMap.h"
#endif

/** Forward declaration of Cmfd class */
class Cmfd;

/**
 * @struct sr_data
 * @brief A sr_data struct represents an SR with a unique SR ID
 *        and a characteristic point that lies within the SR that
 *        can be used to recompute the hierarchical LocalCoords
 *        linked list.
 */
struct sr_data {

  /** The SR ID */
  int _sr_id;

  /** The CMFD Cell */
  int _cmfd_cell;

  /** The Material ID */
  int _mat_id;

  /** Characteristic point in Root Universe that lies in SR */
  Point* _point;

  /** Global numerical centroid in Root Universe */
  Point* _centroid;

  /** Constructor for SR data initializes centroids and points to NULL */
  sr_data() {
    _centroid = NULL;
    _point = NULL;
  }

  /** Destructor for sr_data */
  ~sr_data() {
    if (_point != NULL)
      delete _point;

    if (_centroid != NULL)
      delete _centroid;
  }
};

void reset_auto_ids();


/**
 * @class Geometry Geometry.h "src/Geometry.h"
 * @brief The master class containing references to all geometry-related
 *        objects - Surfaces, Cells, Universes and Lattices - and Materials.
 * @details The primary purpose for the geometry is to serve as a collection
 *          of all geometry-related objects, as well as for ray tracing
 *          of characteristic tracks across the Geometry and computing
 *          SR-to-cell offset maps.
 */
class Geometry {

private:

  /** The boundary conditions at the x-min surface of the bounding box
   *  containing the Geometry. */
  boundaryType _x_min_bc;

  /** The boundary conditions at the y-min surface of the bounding box
   *  containing the Geometry. */
  boundaryType _y_min_bc;

  /** The boundary conditions at the x-max surface of the bounding box
   *  containing the Geometry. */
  boundaryType _x_max_bc;

  /** The boundary conditions at the y-max surface of the bounding box
   *  containing the Geometry. */
  boundaryType _y_max_bc;

  /** An map of SR key hashes to unique sr_data structs */
  ParallelHashMap<std::string, sr_data*> _SR_keys_map;

  /** An vector of SR key hashes indexed by SR ID */
  std::vector<std::string> _SRs_to_keys;

  /* The Universe at the root node in the CSG tree */
  Universe* _root_universe;

  /** A CMFD object pointer */
  Cmfd* _cmfd;

  /* A map of all Material in the Geometry for optimization purposes */
  std::map<int, Material*> _all_materials;

  Cell* findFirstCell(LocalCoords* coords);
  Cell* findNextCell(LocalCoords* coords);

public:

  Geometry();
  virtual ~Geometry();

  /* Get parameters */
  double getWidthX();
  double getWidthY();
  double getWidthZ();
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
  Universe* getRootUniverse();
  int getNumSRs();
  int getNumEnergyGroups();
  int getNumMaterials();
  int getNumCells();
  std::map<int, Material*> getAllMaterials();
  std::map<int, Surface*> getAllSurfaces();
  std::map<int, Cell*> getAllCells();
  std::map<int, Cell*> getAllMaterialCells();
  std::map<int, Universe*> getAllUniverses();
  void setRootUniverse(Universe* root_universe);

  Cmfd* getCmfd();
  std::vector<std::string>& getSRsToKeys();
  int getSRId(LocalCoords* coords);
  Point* getSRPoint(int sr_id);
  Point* getSRCentroid(int sr_id);
  std::string getSRKey(LocalCoords* coords);
  ParallelHashMap<std::string, sr_data*>& getSRKeysMap();

  /* Set parameters */
  void setCmfd(Cmfd* cmfd);
  void setSRCentroid(int sr, Point* centroid);

  /* Find methods */
  Cell* findCellContainingCoords(LocalCoords* coords);
  Material* findSRMaterial(int sr_id);
  int findSRId(LocalCoords* coords);
  Cell* findCellContainingSR(int sr_id);

  /* Other worker methods */
  void subdivideCells();
  void initializeSRs(bool neighbor_cells=false);
  void segmentize(Track* track);
  void initializeSRVectors();
  void computeFissionability(Universe* univ=NULL);

  std::string toString();
  void printString();
  void initializeCmfd();
  bool withinBounds(LocalCoords* coords);
};

#endif /* GEOMETRY_H_ */
