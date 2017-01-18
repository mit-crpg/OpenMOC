#include "../../../src/CPUSolver.h"
#include "../../../src/CPULSSolver.h"
#include "../../../src/log.h"
#include <array>
#include <iostream>

int main(int argc, char* argv[]) {

#ifdef MPIx
  MPI_Init(&argc, &argv);
  log_set_ranks(MPI_COMM_WORLD);
#endif

  /* Define geometry to load */
  std::string file = "beavrs.geo";

  /* Define simulation parameters */
  #ifdef OPENMP
  int num_threads = omp_get_num_procs();
  #else
  int num_threads = 1;
  #endif
  double azim_spacing = 0.2;
  int num_azim = 4;
  double polar_spacing = 0.5;
  int num_polar = 2;
  double tolerance = 1e-6;
  int max_iters = 1400;

  /* Create CMFD lattice */
  Cmfd cmfd;
  //cmfd.setLatticeStructure(23*4, 23*4, 4);
  cmfd.setLatticeStructure(17, 17, 2);
  cmfd.setKNearest(1);

  /* Load the geometry */
  log_printf(NORMAL, "Creating geometry...");
  Geometry geometry;
  geometry.loadFromFile(file);
  geometry.setCmfd(&cmfd);
#ifdef MPIx
  geometry.setDomainDecomposition(2, 2, 2, MPI_COMM_WORLD);
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
  track_generator.setSegmentFormation(OTF_STACKS);
  track_generator.generateTracks();

  /* Run simulation */
  CPULSSolver solver(&track_generator);
  solver.setNumThreads(num_threads);
  solver.setConvergenceThreshold(tolerance);
  solver.computeEigenvalue(max_iters);
  solver.printTimerReport();

  log_printf(TITLE, "Finished");
#ifdef MPIx
  MPI_Finalize();
#endif
  return 0;
}
