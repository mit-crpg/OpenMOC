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
#include "Mesh.h"
#include "Surface.h"
#endif


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

  /** An array of Cell IDs indexed by FSR IDs */
  int* _FSRs_to_cells;

  /** An array of Material UIDs indexed by FSR IDs */
  int* _FSRs_to_material_UIDs;

  /** An array of Material UIDs indexed by FSR IDs */
  int* _FSRs_to_material_IDs;

  /** The maximum Track segment length in the Geometry */
  double _max_seg_length;

  /** The minimum Track segment length in the Geometry */
  double _min_seg_length;

  /* The Universe at the root node in the CSG tree */
  Universe* _root_universe;

  /** A CMFD Mesh object pointer */
  Mesh* _mesh;

  Cell* findFirstCell(LocalCoords* coords, double angle);
  Cell* findNextCell(LocalCoords* coords, double angle);
  Cell* findCell(Universe* univ, int fsr_id);

public:

  Geometry(Mesh* mesh=NULL);
  virtual ~Geometry();

  double getWidth();
  double getHeight();
  double getMinX();
  double getMaxX();
  double getMinY();
  double getMaxY();
  double getMinZ();
  double getMaxZ();
  boundaryType getBCTop();
  boundaryType getBCBottom();
  boundaryType getBCLeft();
  boundaryType getBCRight();
  Universe* getRootUniverse();
  int getNumFSRs();
  int getNumEnergyGroups();
  int getNumMaterials();
  int getNumCells();
  std::map<int, Material*> getAllMaterials();
  int* getFSRtoCellMap();
  int* getFSRtoMaterialMap();
  double getMaxSegmentLength();
  double getMinSegmentLength();
  Mesh* getMesh();

  void setRootUniverse(Universe* root_universe);

  Cell* findCellContainingCoords(LocalCoords* coords);
  CellBasic* findCellContainingFSR(int fsr_id);
  int findFSRId(LocalCoords* coords);
  void subdivideCells();
  void initializeFlatSourceRegions();
  void segmentize(Track* track);
  void computeFissionability(Universe* univ=NULL);

  std::string toString();
  void printString();

  void initializeMesh();
  void findFSRsInCell(Universe* univ, int cell_num, int* fsr_id);
  void defineMesh(Mesh* mesh, Universe* univ, int depth,
                  int* meshCellNum, int row, bool base, int fsr_id);
  int nextLatticeHeight(Universe* univ);
  void findMeshHeight(Universe* univ, int* height, int depth);
  void findMeshWidth(Universe* univ, int* width, int depth);
  int findMeshDepth(Universe* univ, int mesh_level);
};

#endif /* GEOMETRY_H_ */
