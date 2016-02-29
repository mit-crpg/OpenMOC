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
 * @struct fsr_data
 * @brief A fsr_data struct represents an FSR with a unique FSR ID
 *        and a characteristic point that lies within the FSR that
 *        can be used to recompute the hierarchical LocalCoords
 *        linked list.
 */
struct fsr_data {

  /** The FSR ID */
  int _fsr_id;

  /** The CMFD Cell */
  int _cmfd_cell;

  /** The Material ID */
  int _mat_id;

  /** Characteristic point in Root Universe that lies in FSR */
  Point* _point;

  /** Global numerical centroid in Root Universe */
  Point* _centroid;

  /** Constructor for FSR data initializes centroids and points to NULL */
  fsr_data() {
    _centroid = NULL;
    _point = NULL;
  }

  /** Destructor for fsr_data */
  ~fsr_data() {
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
 *          FSR-to-cell offset maps.
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

  /** An map of FSR key hashes to unique fsr_data structs */
  ParallelHashMap<std::string, fsr_data*> _FSR_keys_map;

  /** An vector of FSR key hashes indexed by FSR ID */
  std::vector<std::string> _FSRs_to_keys;

  /* The Universe at the root node in the CSG tree */
  Universe* _root_universe;

  /** A CMFD object pointer */
  Cmfd* _cmfd;

  Point** _centroids;

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
  int getNumFSRs();
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
  std::vector<std::string> & getFSRsToKeys();
  int getFSRId(LocalCoords* coords);
  Point* getFSRPoint(int fsr_id);
  Point* getFSRCentroid(int fsr_id);
  std::string getFSRKey(LocalCoords* coords);
  ParallelHashMap<std::string, fsr_data*> & getFSRKeysMap();

  /* Set parameters */
  void setCmfd(Cmfd* cmfd);
  void setFSRCentroid(int fsr, Point* centroid);

  /* Find methods */
  Cell* findCellContainingCoords(LocalCoords* coords);
  Material* findFSRMaterial(int fsr_id);
  int findFSRId(LocalCoords* coords);
  Cell* findCellContainingFSR(int fsr_id);

  /* Other worker methods */
  void subdivideCells();
  void initializeFSRs(bool neighbor_cells=false);
  void segmentize(Track* track);
  void initializeFSRVectors();
  void computeFissionability(Universe* univ=NULL);

  std::string toString();
  void printString();
  void initializeCmfd();
  bool withinBounds(LocalCoords* coords);
};

#endif /* GEOMETRY_H_ */
