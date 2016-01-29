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
#include "Progress.h"
#include <limits>
#include <sys/types.h>
#include <sys/stat.h>
#include <sstream>
#include <string>
#include <set>
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

  /** Constructor for FSR data object */
  fsr_data() : _fsr_id(0), _cmfd_cell(0), _mat_id(0), _point(NULL),
    _centroid(NULL){}

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

  /** Destructor for fsr_data */
  ~fsr_data() {
    if (_point != NULL)
      delete _point;

    if (_centroid != NULL)
      delete _centroid;
  }
};


/**
 * @struct ExtrudedFSR
 * @brief An ExtrudedFSR struct represents a FSR region in the superposition
 *        plane for axial on-the-fly ray tracing. It contains a characteristic
 *        point that lies within the FSR, an axial mesh, and an array of 3D
 *        FSR IDs contained within the extruded region along with their
 *        corresponding materials.
 */
struct ExtrudedFSR {

  /** Constructor for ExtrudedFSR object */
  ExtrudedFSR() : _mesh(NULL), _fsr_id(0), _fsr_ids(NULL), _materials(NULL),
    _num_fsrs(0), _coords(NULL){}

  /** Array defining the axial mesh */
  FP_PRECISION* _mesh;

  /** Axial extruded FSR ID */
  int _fsr_id;

  /** Array of 3D FSR IDs */
  int* _fsr_ids;

  /** Array of material pointers for each FSR */
  Material** _materials;

  /** Number of FSRs in the axially extruded FSR */
  size_t _num_fsrs;

  /** Coordinates inside the FSR */
  LocalCoords* _coords;
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

  bool _solve_3D;

  /** The boundary conditions at the top of the bounding box containing
   *  the Geometry. False is for vacuum and true is for reflective BCs. */
  boundaryType _x_min_bc;

  /** The boundary conditions at the top of the bounding box containing
   *  the Geometry. False is for vacuum and true is for reflective BCs. */
  boundaryType _x_max_bc;

  /** The boundary conditions at the top of the bounding box containing
   *  the Geometry. False is for vacuum and true is for reflective BCs. */
  boundaryType _y_min_bc;

  /** The boundary conditions at the top of the bounding box containing
   *  the Geometry. False is for vacuum and true is for reflective BCs. */
  boundaryType _y_max_bc;

  /** The boundary conditions at the top of the bounding box containing
   *  the Geometry. False is for vacuum and true is for reflective BCs. */
  boundaryType _z_min_bc;

  /** The boundary conditions at the top of the bounding box containing
   *  the Geometry. False is for vacuum and true is for reflective BCs. */
  boundaryType _z_max_bc;

  /** An map of FSR key hashes to unique fsr_data structs */
  ParallelHashMap<std::size_t, fsr_data*> _FSR_keys_map;
  ParallelHashMap<std::size_t, ExtrudedFSR*> _extruded_FSR_keys_map;

  /** An vector of FSR key hashes indexed by FSR ID */
  std::vector<std::size_t> _FSRs_to_keys;

  /** A vector of Material IDs indexed by FSR IDs */
  std::vector<int> _FSRs_to_material_IDs;

  /** A vector of ExtrudedFSR pointers indexed by extruded FSR ID */
  std::vector<ExtrudedFSR*> _extruded_FSR_lookup;

  /* The Universe at the root node in the CSG tree */
  Universe* _root_universe;

  /** A CMFD object pointer */
  Cmfd* _cmfd;

  /* A map of all Material in the Geometry for optimization purposes */
  std::map<int, Material*> _all_materials;

  Cell* findFirstCell(LocalCoords* coords, double azim, double polar=M_PI_2);
  Cell* findNextCell(LocalCoords* coords, double azim, double polar=M_PI_2);

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
  boundaryType getMinZBoundaryType();
  boundaryType getMaxZBoundaryType();
  Universe* getRootUniverse();
  int getNumFSRs();
  int getNumEnergyGroups();
  int getNumMaterials();
  int getNumCells();
  std::map<int, Material*> getAllMaterials();
  std::map<int, Cell*> getAllMaterialCells();
  std::vector<FP_PRECISION> getUniqueZHeights();
  std::vector<FP_PRECISION> getUniqueZPlanes();
  void setRootUniverse(Universe* root_universe);

  Cmfd* getCmfd();
  std::vector<std::size_t>* getFSRsToKeys();
  std::vector<int>* getFSRsToMaterialIDs();
  int getFSRId(LocalCoords* coords);
  Point* getFSRPoint(int fsr_id);
  Point* getFSRCentroid(int fsr_id);
  int getCmfdCell(int fsr_id);
  ExtrudedFSR* getExtrudedFSR(int extruded_fsr_id);
  std::string getFSRKey(LocalCoords* coords);
  ParallelHashMap<std::size_t, fsr_data*>* getFSRKeysMap();

  /* Set parameters */
  void setFSRsToMaterialIDs(std::vector<int>* FSRs_to_material_IDs);
  void setFSRsToKeys(std::vector<std::size_t>* FSRs_to_keys);
  void setCmfd(Cmfd* cmfd);
  void setFSRCentroid(int fsr, Point* centroid);
  void setFSRKeysMap(ParallelHashMap<std::size_t, fsr_data*>* FSR_keys_map);

  /* Find methods */
  Cell* findCellContainingCoords(LocalCoords* coords);
  Material* findFSRMaterial(int fsr_id);
  int findFSRId(LocalCoords* coords);
  int findExtrudedFSR(LocalCoords* coords);
  Cell* findCellContainingFSR(int fsr_id);

  /* Other worker methods */
  void subdivideCells();
  void initializeAxialFSRs(std::vector<double> global_z_mesh);
  void initializeFlatSourceRegions();
  void segmentize2D(Track2D* track, double z_coord);
  void segmentize3D(Track3D* track);
  void segmentizeExtruded(Track* flattened_track,
      std::vector<double> z_coords);
  void initializeFSRVectors();
  void computeFissionability(Universe* univ=NULL);

  std::string toString();
  void printString();
  void initializeCmfd();
  bool withinBounds(LocalCoords* coords);
};

#endif /* GEOMETRY_H_ */
