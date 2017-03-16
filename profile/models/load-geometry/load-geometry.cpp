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
  //std::string file = "few-pins-reduced-ref-v2.geo";
  //std::string file = "full-assembly-final.geo";
  //std::string file = "beavrs-2D-v2-NULL-corr.geo";
  //std::string file = "beavrs-2D-v2.geo";
  //std::string file = "beavrs-2D-v3.geo";
  std::string file = "full-core-beavrs-new.geo";

  /* Define simulation parameters */
  #ifdef OPENMP
  int num_threads = omp_get_num_procs();
  #else
  int num_threads = 1;
  #endif
 
  double azim_spacing = 0.2;
  int num_azim = 8; //FIXME 32
  double polar_spacing = 1.0;
  int num_polar = 6;

  /*
  double azim_spacing = 0.1;
  int num_azim = 32;
  double polar_spacing = 5.0; // 1.0
  int num_polar = 6;
  */

  double tolerance = 1e-7;
  int max_iters = 250;

  /* Create CMFD lattice */
  Cmfd cmfd;
  cmfd.useAxialInterpolation(true);
  cmfd.setLatticeStructure(17*17, 17*17, 200);
  cmfd.setKNearest(1);
  std::vector<std::vector<int> > cmfd_group_structure =
      get_group_structure(70, 8);
  cmfd.setGroupStructure(cmfd_group_structure);
  cmfd.setCMFDRelaxationFactor(0.7);
  cmfd.setSORRelaxationFactor(1.6);
  cmfd.useFluxLimiting(true);

  /* Load the geometry */
  log_printf(NORMAL, "Creating geometry...");
  Geometry geometry;
  geometry.loadFromFile(file);
  //geometry.manipulateXS();

  //FIXME geometry.setOverlaidMesh(2.0);
  geometry.setCmfd(&cmfd);
  log_printf(NORMAL, "Pitch = %8.6e", geometry.getMaxX() - geometry.getMinX());
  log_printf(NORMAL, "Height = %8.6e", geometry.getMaxZ() - geometry.getMinZ());
#ifdef MPIx
  geometry.setDomainDecomposition(17, 17, 4, MPI_COMM_WORLD);
  //geometry.setNumDomainModules(17,17,1);
  //geometry.setDomainDecomposition(17, 17, 1, MPI_COMM_WORLD); // FIXME 23
#else
  //geometry.setNumDomainModules(2,2,4);
#endif
  geometry.initializeFlatSourceRegions();

  /* Create the track generator */
  log_printf(NORMAL, "Initializing the track generator...");
  TrackGenerator3D track_generator(&geometry, num_azim, num_polar, azim_spacing,
                                   polar_spacing);
  track_generator.setTrackGenerationMethod(MODULAR_RAY_TRACING);
  track_generator.setNumThreads(num_threads);
  track_generator.setSegmentFormation(OTF_STACKS);
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
  solver.setVerboseIterationReport();
  solver.setConvergenceThreshold(tolerance);
  //solver.correctXS();
  solver.setCheckXSLogLevel(INFO);
  solver.computeEigenvalue(max_iters);
  solver.printTimerReport();

  Lattice mesh_lattice;
  Mesh mesh(&solver);
  mesh.createLattice(17*17, 17*17, 1);
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
