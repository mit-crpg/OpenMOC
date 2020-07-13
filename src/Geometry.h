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
#include <fstream>
#include <string>
#include <set>
#include <omp.h>
#include <functional>
#include "ParallelHashMap.h"
#endif

#ifdef MPIx
#include <mpi.h>
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
  long _fsr_id;

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
  double* _mesh;

  /** Axial extruded FSR ID */
  int _fsr_id;

  /** Array of 3D FSR IDs */
  long* _fsr_ids;

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

  /** An map of FSR key hashes to unique fsr_data structs */
  ParallelHashMap<std::string, fsr_data*> _FSR_keys_map;
  ParallelHashMap<std::string, ExtrudedFSR*> _extruded_FSR_keys_map;

  /** An vector of FSR key hashes indexed by FSR ID */
  std::vector<std::string> _FSRs_to_keys;

  /** An vector of FSR centroids indexed by FSR ID */
  std::vector<Point*> _FSRs_to_centroids;

  /** A boolean indicating whether any centroids have been set */
  bool _contains_FSR_centroids;

  /** A vector of Material IDs indexed by FSR IDs */
  std::vector<int> _FSRs_to_material_IDs;

  /** A vector of ExtrudedFSR pointers indexed by extruded FSR ID */
  std::vector<ExtrudedFSR*> _extruded_FSR_lookup;

  /** An vector of CMFD cell IDs indexed by FSR ID */
  std::vector<int> _FSRs_to_CMFD_cells;

  /* The Universe at the root node in the CSG tree */
  Universe* _root_universe;

  /** A CMFD object pointer */
  Cmfd* _cmfd;

  /** An optional axial mesh overlaid on the Geometry */
  Lattice* _overlaid_mesh;

  /* A map of all Material in the Geometry for optimization purposes */
  std::map<int, Material*> _all_materials;

  /* A vector containing allocated strings for key generation */
  std::vector<std::string> _fsr_keys;

  /* A boolean to know whether geometry is domain decomposed or not */
  bool _domain_decomposed;

  /* A boolean to know whether FSRs were counted */
  bool _domain_FSRs_counted;

  /* Number of domains in the X, Y and Z directions */
  int _num_domains_x;
  int _num_domains_y;
  int _num_domains_z;

  /* Index of the domain in the whole geometry in the X, Y and Z directions */
  int _domain_index_x;
  int _domain_index_y;
  int _domain_index_z;

  /* Lattice object of size 1 that contains the local domain */
  Lattice* _domain_bounds;

  /* Number of FSRs in each domain */
  std::vector<long> _num_domain_FSRs;

#ifdef MPIx
  /* MPI communicator, used to transfer information across domain for both
   * ray-tracing and solving the MOC equations */
  MPI_Comm _MPI_cart;
#endif

  /* Number of modular track laydowns within a domain, in the X, Y and
     Z directions */
  int _num_modules_x;
  int _num_modules_y;
  int _num_modules_z;

  /* Symmetries in X,Y and Z axis used to restrict the computation domain */
  std::vector<bool> _symmetries;

  /* Wether to read bites backwards or forward (for BGQ) */
  bool _twiddle;

  /* Whether the geometry was loaded from a .geo file or created in the input
   * file. This matters for memory de-allocation purposes. */
  bool _loaded_from_file;

  /* Function to find the cell containing the coordinates */
  Cell* findFirstCell(LocalCoords* coords, double azim, double polar=M_PI_2);

public:

  Geometry();
  virtual ~Geometry();

  /* Handle number of ray-tracing modules in a domain */
  void setNumDomainModules(int num_x, int num_y, int num_z);
  int getNumXModules();
  int getNumYModules();
  int getNumZModules();

  /* Handle symmetry axis used to restrict the computation domain */
  void useSymmetry(bool X_symmetry, bool Y_symmetry, bool Z_symmetry);
  bool getSymmetry(int axis);

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
  long getNumFSRs();
  long getNumTotalFSRs();
  int getNumEnergyGroups();
  int getNumMaterials();
  int getNumCells();
  std::map<int, Surface*> getAllSurfaces();
  std::map<int, Material*> getAllMaterials();
  std::map<int, Cell*> getAllCells();
  std::map<int, Cell*> getAllMaterialCells();
  std::map<int, Universe*> getAllUniverses();
  std::vector<double> getUniqueZHeights(bool include_overlaid_mesh = false);
  std::vector<double> getUniqueZPlanes();
  bool isDomainDecomposed();
  bool isRootDomain();
  void getDomainIndexes(int* indexes);
  void getDomainStructure(int* structure);

  /* Assign root universe to geometry */
  void setRootUniverse(Universe* root_universe);

  /* Set up domain decomposition */
#ifdef MPIx
  void setDomainDecomposition(int nx, int ny, int nz, MPI_Comm comm);
  MPI_Comm getMPICart();
#endif

  /* More complex getter methods */
  Cmfd* getCmfd();
  std::vector<std::string>& getFSRsToKeys();
  std::vector<int>& getFSRsToMaterialIDs();
  std::vector<Point*>& getFSRsToCentroids();
  std::vector<int>& getFSRsToCMFDCells();
  std::vector<ExtrudedFSR*>& getExtrudedFSRLookup();
  long getFSRId(LocalCoords* coords, bool err_check=true);
  long getGlobalFSRId(LocalCoords* coords, bool err_check=true);
  Point* getFSRPoint(long fsr_id);
  Point* getFSRCentroid(long fsr_id);
  bool containsFSRCentroids();
  int getCmfdCell(long fsr_id);
  ExtrudedFSR* getExtrudedFSR(int extruded_fsr_id);
  std::string getFSRKey(LocalCoords* coords);
  void getFSRKeyFast(LocalCoords* coords, std::string& key);
  void printToString(std::string& str, int& index, int value);
  int getNumDigits(int number);
  ParallelHashMap<std::string, fsr_data*>& getFSRKeysMap();
  ParallelHashMap<std::string, ExtrudedFSR*>& getExtrudedFSRKeysMap();
#ifdef MPIx
  int getNeighborDomain(int offset_x, int offset_y, int offset_z);
#endif

  /* Setter methods */
  void setCmfd(Cmfd* cmfd);
  void setFSRCentroid(long fsr, Point* centroid);
  void resetContainsFSRCentroids();
  void setOverlaidMesh(double axial_mesh_height, int num_x=0,
                       int num_y=0, int num_radial_domains=0,
                       int* radial_domains=NULL);
  void clearBoundaries();

  /* Find methods */
  Cell* findCellContainingCoords(LocalCoords* coords);
  Material* findFSRMaterial(long fsr_id);
  long findFSRId(LocalCoords* coords);
  int findExtrudedFSR(LocalCoords* coords);
  Cell* findCellContainingFSR(long fsr_id);
  Cell* findNextCell(LocalCoords* coords, double azim, double polar=M_PI_2);

  /* Other worker methods */
  void setNumThreads(int num_threads);
  void reserveKeyStrings(int num_threads);
  void subdivideCells();
  void initializeAxialFSRs(std::vector<double> global_z_mesh);
  void reorderFSRIDs();
  void initializeFlatSourceRegions();
  void segmentize2D(Track* track, double z_coord);
  void segmentize3D(Track3D* track, bool setup=false);
  void segmentizeExtruded(Track* flattened_track,
                          std::vector<double> z_coords);
  void fixFSRMaps();
  void initializeFSRVectors();
  void computeFissionability(Universe* univ=NULL);
  void manipulateXS();
  void loadSPHFactors(double* sph_factors, int num_domains_groups,
                      double* sph_to_domain_ids, int num_sph_domains,
                      const char* domain_type);

  /* Obtain or print information about the geometry */
  std::vector<long> getSpatialDataOnGrid(std::vector<double> dim1,
                                         std::vector<double> dim2,
                                         double offset,
                                         const char* plane,
                                         const char* domain_type);
  std::string toString();
  void printString();
  void printFSRsToFile(const char* plane="xy", int gridsize=1000,
                       double offset=0.0, double* bounds_x = NULL,
                       double* bounds_y = NULL, double* bounds_z = NULL);

  void initializeCmfd();
  void initializeSpectrumCalculator(Cmfd* spectrum_calculator);
  bool withinBounds(LocalCoords* coords);
  bool withinGlobalBounds(LocalCoords* coords);
#ifdef MPIx
  void countDomainFSRs();
  void getLocalFSRId(long global_fsr_id, long &local_fsr_id, int &domain);
#endif
  std::vector<double> getGlobalFSRCentroidData(long global_fsr_id);
  int getDomainByCoords(LocalCoords* coords);
  std::map<Cell*, std::vector<long> > getCellsToFSRs();

  /* Input/output of geometries from/to .geo files */
  void dumpToFile(std::string filename);
  void loadFromFile(std::string filename, bool twiddle=false);
  size_t twiddleRead(int* ptr, size_t size, size_t nmemb, FILE* stream);
  size_t twiddleRead(bool* ptr, size_t size, size_t nmemb, FILE* stream);
  size_t twiddleRead(char* ptr, size_t size, size_t nmemb, FILE* stream);
  size_t twiddleRead(universeType* ptr, size_t size, size_t nmemb, FILE* stream);
  size_t twiddleRead(cellType* ptr, size_t size, size_t nmemb, FILE* stream);
  size_t twiddleRead(surfaceType* ptr, size_t size, size_t nmemb, FILE* stream);
  size_t twiddleRead(boundaryType* ptr, size_t size, size_t nmemb, FILE* stream);
  size_t twiddleRead(double* ptr, size_t size, size_t nmemb, FILE* stream);
  size_t twiddleRead(long* ptr, size_t size, size_t nmemb, FILE* stream);
};

#endif /* GEOMETRY_H_ */
