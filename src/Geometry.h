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

  /** A CMFD Mesh object pointer */
  Mesh* _mesh;

  void initializeCellFillPointers();

  Cell* findFirstCell(LocalCoords* coords, double angle);
  Cell* findNextCell(LocalCoords* coords, double angle);
  Cell* findCell(Universe* univ, int fsr_id);

public:

  Geometry(Mesh* mesh=NULL);
  virtual ~Geometry();

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
  int* getFSRtoCellMap();
  int* getFSRtoMaterialMap();
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
  Mesh* getMesh();

  void addMaterial(Material* material);
  void addSurface(Surface* surface);
  void addCell(Cell *cell);
  void addUniverse(Universe* universe);
  void addLattice(Lattice* lattice);

  void removeMaterial(int id);
  void removeCell(int id);
  void removeUniverse(int id);
  void removeLattice(int id);

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
