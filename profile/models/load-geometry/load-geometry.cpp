#include "../../../src/CPUSolver.h"
#include "../../../src/CPULSSolver.h"
#include "../../../src/log.h"
#include "../../../src/Mesh.h"
#include <array>
#include <iostream>
#include "helper-code/group-structures.h"

int main(int argc,  char* argv[]) {

#ifdef MPIx
  MPI_Init(&argc, &argv);
  log_set_ranks(MPI_COMM_WORLD);
#endif

  /* Define geometry to load */
  //std::string file = "full-axial-simple-assembly.geo";
  //std::string file = "full-axial-w-refs.geo";
  //std::string file = "assembly-1.6-70grp-discr.geo";
  //std::string file = "full-axial-detail.geo";
  //std::string file = "full-detail-pin-cell.geo";
  //std::string file = "pin-cell-full-axial-no-refs.geo";
  std::string file = "pin-cell-simple.geo";

  /* Define simulation parameters */
  #ifdef OPENMP
  int num_threads = omp_get_num_procs();
  #else
  int num_threads = 1;
  #endif
  
  double azim_spacing = 0.5;
  int num_azim = 8;
  double polar_spacing = 1.; // 1.0
  int num_polar = 2;
  double tolerance = 1e-4;
  int max_iters = 50;

  /* Create CMFD lattice */
  int num_fsrs = 200;
  Cmfd cmfd;
  cmfd.setLatticeStructure(1, 1, 20);
  cmfd.useAxialInterpolation(true);
  //cmfd.setLatticeStructure(17, 17, 230);
  cmfd.setKNearest(3);
  std::vector<std::vector<int> > cmfd_group_structure = 
      get_group_structure(2, 2);
  cmfd.setGroupStructure(cmfd_group_structure);
  cmfd.setCMFDRelaxationFactor(1.0);
  cmfd.setSORRelaxationFactor(1.0);
  //cmfd.setCMFDRelaxationFactor(0.5);

  /* Load the geometry */
  log_printf(NORMAL, "Creating geometry...");
  Geometry geometry2;
  geometry2.loadFromFile(file);

  XPlane xmin(17*geometry2.getMinX());
  XPlane xmax(17*geometry2.getMaxX());
  YPlane ymin(17*geometry2.getMinY());
  YPlane ymax(17*geometry2.getMaxY());
  ZPlane zmin(geometry2.getMinZ());
  ZPlane zmax(geometry2.getMaxZ());

  xmin.setBoundaryType(REFLECTIVE);
  ymin.setBoundaryType(REFLECTIVE);
  zmin.setBoundaryType(VACUUM);
  xmax.setBoundaryType(REFLECTIVE);
  ymax.setBoundaryType(REFLECTIVE);
  zmax.setBoundaryType(VACUUM);
  

  Cell* fuel = new Cell();
  Material* mat = new Material(10000, "fuel");
  fuel->setFill(mat);
  fuel->addSurface(+1, &xmin);
  fuel->addSurface(-1, &xmax);
  fuel->addSurface(+1, &ymin);
  fuel->addSurface(-1, &ymax);
  fuel->addSurface(+1, &zmin);
  fuel->addSurface(-1, &zmax);
  
  Universe* root_universe = new Universe();
  root_universe->addCell(fuel);
  Geometry geometry;
  geometry.setRootUniverse(root_universe);

  double sigma_t[2] = {0.25, 1.0};
  double sigma_s[4] = {0.235, 0.015, 0.0, 0.9};
  double nu_sigma_f[2] = {0.0, 0.15};
  double sigma_f[2] = {0.0, 0.15 / 2.4};
  double chi[2] = {1.0, 0.0};

  std::map<int, Material*> materials = geometry.getAllMaterials();
  int ng = 2;
  std::map<int, Material*>::iterator material_iter;
  for (material_iter = materials.begin();
        material_iter != materials.end(); ++material_iter) {
    int key = material_iter->first;
    mat->setNumEnergyGroups(ng);
    for (int i=0; i < ng; i++) {
        mat->setSigmaTByGroup(sigma_t[i], i+1);
        mat->setSigmaFByGroup(sigma_f[i], i+1);
        mat->setNuSigmaFByGroup(nu_sigma_f[i], i+1);
        mat->setChiByGroup(chi[i], i+1);
        for (int j=0; j < ng; j++)
          mat->setSigmaSByGroup(sigma_s[2*i+j], i+1, j+1);
      }
  }


  std::map<int, Cell*> cells = geometry.getRootUniverse()->getCells();
  std::map<int, Cell*>::iterator c_iter;
  std::map<int, surface_halfspace*>::iterator s_iter;
  Surface* surf;
  int halfspace;
  double z_min = geometry.getMinZ();
  double z_max = geometry.getMaxZ();
  for (c_iter = cells.begin(); c_iter != cells.end(); ++c_iter) {
    std::map<int, surface_halfspace*> surfs = c_iter->second->getSurfaces();

    for (s_iter = surfs.begin(); s_iter != surfs.end(); ++s_iter) {
      surf = s_iter->second->_surface;
      halfspace = s_iter->second->_halfspace;

      /*
      if (surf->getSurfaceType() == ZPLANE) {
        if (surf->getBoundaryType() == VACUUM) {
            surf->setBoundaryType(REFLECTIVE);
        }
      }
      */
    }
  }
  geometry.getRootUniverse()->calculateBoundaries();


  geometry.setAxialMesh(400.0/num_fsrs);
  geometry.setCmfd(&cmfd);
  log_printf(NORMAL, "Pitch = %8.6e", geometry.getMaxX() - geometry.getMinX());
#ifdef MPIx
  ///geometry.setDomainDecomposition(2, 2, 23, MPI_COMM_WORLD); // FIXME 23
  //geometry.setNumDomainModules(2,2,2);
#else
  geometry.setNumDomainModules(2,2,2);
#endif
  geometry.setNumDomainModules(1,1,num_fsrs);
  geometry.initializeFlatSourceRegions();

  /* Create the track generator */
  log_printf(NORMAL, "Initializing the track generator...");
  TrackGenerator3D track_generator(&geometry, num_azim, num_polar, azim_spacing,
                                   polar_spacing);
  track_generator.setTrackGenerationMethod(MODULAR_RAY_TRACING);
  track_generator.setNumThreads(num_threads);
  track_generator.setSegmentFormation(OTF_TRACKS);
  /*
  std::vector<FP_PRECISION> zones;
  zones.push_back(0.0);
  zones.push_back(20.0);
  zones.push_back(380.0);
  zones.push_back(400.0);
  track_generator.setSegmentationZones(zones);
  */
  track_generator.generateTracks();

  /* Run simulation */
  CPUSolver solver(&track_generator);
  solver.setNumThreads(num_threads);
  solver.setConvergenceThreshold(tolerance);
  solver.correctXS();
  solver.setCheckXSLogLevel(WARNING);
  solver.computeEigenvalue(max_iters);
  solver.printTimerReport();

  Lattice mesh_lattice;
  Mesh mesh(&solver);
  //FIXME mesh.createLattice(17, 17, 230);
  mesh.createLattice(1, 1, num_fsrs);
  Vector3D rx_rates = mesh.getFormattedReactionRates(FISSION_RX);

  int my_rank = 0;
#ifdef MPIx
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#endif
  if (my_rank == 0) {
    for (int k=0; k < rx_rates.at(0).at(0).size(); k++) {
      std::cout << " -------- z = " << k << " ----------" << std::endl;
      for (int j=0; j < rx_rates.at(0).size(); j++) {
        for (int i=0; i < rx_rates.size(); i++) {
          std::cout << rx_rates.at(i).at(j).at(k) << " ";
        }
        std::cout << std::endl;
      }
    }
  }

  log_printf(TITLE, "Finished");
#ifdef MPIx
  MPI_Finalize();
#endif
  return 0;
}
