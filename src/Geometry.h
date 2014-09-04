/**
 * @file Geometry.h
 * @brief The Geometry class.
 * @date January 19, 2012
 * @author William Boyd, MIT, Course 22 (wboyd@mit.edu)
 */

#ifndef GEOMETRY_H_
#define GEOMETRY_H_

#ifdef __cplusplus
#include <limits.h>
#include <limits>
#include <sys/types.h>
#include <sys/stat.h>
#include "LocalCoords.h"
#include "Track.h"
#include "Surface.h"
#include "Cmfd.h"
#include <sstream>
#include <string>
#include <omp.h>
#endif

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

  /** Characteristic point in Universe 0 that lies in FSR */
  Point* _point;

};

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

  omp_lock_t* _num_FSRs_lock;

  /** The minimum point along the x-axis contained by Geometry in cm */
  double _x_min;

  /** The minimum point along the y-axis contained by Geometry in cm */
  double _y_min;

  /** The maximum point along the x-axis contained by Geometry in cm */
  double _x_max;

  /** The maximum point along the y-axis contained by Geometry in cm */
  double _y_max;

  /** The boundary conditions at the top of the bounding box containing
   *  the Geometry. False is for vacuum and true is for reflective BCs. */
  boundaryType _top_bc;

  /** The boundary conditions at the top of the bounding box containing
   *  the Geometry. False is for vacuum and true is for reflective BCs. */
  boundaryType _bottom_bc;

  /** The boundary conditions at the top of the bounding box containing
   *  the Geometry. False is for vacuum and true is for reflective BCs. */
  boundaryType _left_bc;

  /** The boundary conditions at the top of the bounding box containing
   *  the Geometry. False is for vacuum and true is for reflective BCs. */
  boundaryType _right_bc;

  /** The total number of FSRs in the Geometry */
  int _num_FSRs;

  /** The number of energy groups for each Material's nuclear data */
  int _num_groups;

  /** An map of FSR key hashes to unique fsr_data structs */
  std::map<std::size_t, fsr_data> _FSR_keys_map;

  /** An vector of FSR key hashes indexed by FSR ID */
  std::vector<std::size_t> _FSRs_to_keys;

  /** A vector of Material IDs indexed by FSR IDs */
  std::vector<int> _FSRs_to_material_IDs;

  /** The maximum Track segment length in the Geometry */
  double _max_seg_length;

  /** The minimum Track segment length in the Geometry */
  double _min_seg_length;

  /** A std::map of Material IDs (keys) to Material object pointers (values) */
  std::map<int, Material*> _materials;

  /** A std::map of Surface IDs (keys) to Surface object pointers (values) */
  std::map<int, Surface*> _surfaces;

  /** A std::map of Cell IDs (keys) to Cell object pointers (values) */
  std::map<int, Cell*> _cells;

  /** A std::map of Universe IDs (keys) to Universe object pointers (values) */
  std::map<int, Universe*> _universes;

  /** A std::map of Lattice IDs (keys) to Lattice object pointers (values) */
  std::map<int, Lattice*> _lattices;

  /** A CMFD object pointer */
  Cmfd* _cmfd;

  void initializeCellFillPointers();  
  CellBasic* findFirstCell(LocalCoords* coords, double angle);
  CellBasic* findNextCell(LocalCoords* coords, double angle);


public:

  Geometry(Cmfd* cmfd=NULL);
  virtual ~Geometry();

  /* get parameters */
  double getWidth();
  double getHeight();
  double getXMin();
  double getXMax();
  double getYMin();
  double getYMax();
  boundaryType getBCTop();
  boundaryType getBCBottom();
  boundaryType getBCLeft();
  boundaryType getBCRight();
  int getNumFSRs();
  int getNumEnergyGroups();
  int getNumMaterials();
  int getNumCells();

  double getMaxSegmentLength();
  double getMinSegmentLength();
  std::map<int, Material*> getMaterials();
  Material* getMaterial(int id);
  Surface* getSurface(int id);
  Cell* getCell(int id);
  CellBasic* getCellBasic(int id);
  CellFill* getCellFill(int id);
  Universe* getUniverse(int id);
  Lattice* getLattice(int id);
  Cmfd* getCmfd();
  std::map<std::size_t, fsr_data> getFSRKeysMap();
  std::vector<std::size_t> getFSRsToKeys();
  std::vector<int> getFSRsToMaterials();
  int getFSRId(LocalCoords* coords);
  Point* getFSRPoint(int fsr_id);
  std::string getFSRKey(LocalCoords* coords);

  /* set parameters */
  void setFSRKeysMap(std::map<std::size_t, fsr_data> FSR_keys_map);
  void setFSRsToMaterials(std::vector<int> FSRs_to_material_IDs);
  void setFSRsToKeys(std::vector<std::size_t> FSRs_to_keys);
  void setNumFSRs(int num_fsrs);

  /* add object methods */
  void addMaterial(Material* material);
  void addSurface(Surface* surface);
  void addCell(Cell *cell);
  void addUniverse(Universe* universe);
  void addLattice(Lattice* lattice);

  /* remove object methods */
  void removeMaterial(int id);
  void removeCell(int id);
  void removeUniverse(int id);
  void removeLattice(int id);

  /* find methods */
  CellBasic* findCellContainingCoords(LocalCoords* coords);
  Material* findMaterialContainingFSR(int fsr_id);
  int findFSRId(LocalCoords* coords);

  /* Other worker methods */
  void subdivideCells();
  void initializeFlatSourceRegions();
  void segmentize(Track* track);
  void computeFissionability(Universe* univ=NULL);
  std::string toString();
  void printString();
  void initializeCmfd();

};

#endif /* GEOMETRY_H_ */
