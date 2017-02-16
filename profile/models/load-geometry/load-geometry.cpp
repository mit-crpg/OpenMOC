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
  std::string file = "full-axial-detail.geo";

  /* Define simulation parameters */
  #ifdef OPENMP
  int num_threads = omp_get_num_procs();
  #else
  int num_threads = 1;
  #endif
  
  double azim_spacing = 0.05;
  int num_azim = 64;
  double polar_spacing = 0.25;
  int num_polar = 14;
  double tolerance = 1e-4;
  int max_iters = 50;

  /* Create CMFD lattice */
  Cmfd cmfd;
  //cmfd.setLatticeStructure(17, 17, 200);
  //cmfd.useAxialInterpolation(true);
  cmfd.setLatticeStructure(17, 17, 230);
  cmfd.setKNearest(1);
  std::vector<std::vector<int> > cmfd_group_structure = 
      get_group_structure(70,8);
  cmfd.setGroupStructure(cmfd_group_structure);
  cmfd.setCMFDRelaxationFactor(0.8);
  //cmfd.setCMFDRelaxationFactor(0.5);

  /* Load the geometry */
  log_printf(NORMAL, "Creating geometry...");
  Geometry geometry;
  geometry.loadFromFile(file);
  geometry.setAxialMesh(0.5);
  
  //FIXME
  /*
  std::map<int, Material*> materials = geometry.getAllMaterials();
  int ng = geometry.getNumEnergyGroups();
  //FIXME groups 37 and 39 are the issue in helium
  std::string file2 = "pin-cell.geo"; // (discr)
  Geometry geometry2;
  geometry2.loadFromFile(file2);
  std::map<int, Material*> materials2 = geometry2.getAllMaterials();
  std::map<int, Material*>::iterator material_iter;
  for (material_iter = materials.begin();
      material_iter != materials.end(); ++material_iter) {
    int key = material_iter->first;
    Material* mat = materials[key];
    Material* mat2;
    if (materials2.find(key) == materials2.end())
        mat2 = materials[10014];
    else
        mat2 = materials2[key];
    for (int i=0; i < ng; i++) {
        mat->setSigmaTByGroup(mat2->getSigmaTByGroup(i+1), i+1);
        mat->setSigmaFByGroup(mat2->getSigmaFByGroup(i+1), i+1);
        mat->setNuSigmaFByGroup(mat2->getNuSigmaFByGroup(i+1), i+1);
        mat->setChiByGroup(mat2->getChiByGroup(i+1), i+1);
        for (int j=0; j < ng; j++)
          mat->setSigmaSByGroup(mat2->getSigmaSByGroup(i+1,j+1), i+1, j+1);
      }
  }
  */

  geometry.setCmfd(&cmfd);
#ifdef MPIx
  geometry.setDomainDecomposition(2, 2, 23, MPI_COMM_WORLD); // FIXME 23
  //geometry.setNumDomainModules(2,2,2);
#else
  geometry.setNumDomainModules(2,2,2);
#endif
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
  CPULSSolver solver(&track_generator);
  //solver.useExponentialIntrinsic();
  solver.setNumThreads(num_threads);
  solver.setConvergenceThreshold(tolerance);
  solver.correctXS();
  solver.setCheckXSLogLevel(WARNING);
  solver.computeEigenvalue(max_iters);
  solver.printTimerReport();

  Lattice mesh_lattice;
  Mesh mesh(&solver);
  mesh.createLattice(17, 17, 230);
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
