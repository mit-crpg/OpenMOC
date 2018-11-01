#include "../../../src/CPUSolver.h"
#include "../../../src/CPULSSolver.h"
#include "../../../src/log.h"
#include "../../../src/Mesh.h"
#include <array>
#include <iostream>
#include <fenv.h>

int main(int argc,  char* argv[]) {

#ifdef MPIx
  MPI_Init(&argc, &argv);
  log_set_ranks(MPI_COMM_WORLD);
#endif

  set_line_length(120);
  std::string file = "non-uniform-lattice.geo";
  
  double start_time = omp_get_wtime();

  /* Define simulation parameters */
  #ifdef OPENMP
  int num_threads = 1;
  //omp_get_num_procs();
  #else
  int num_threads = 1;
  #endif
 
  double azim_spacing = 0.05;
  int num_azim = 64; 
  double polar_spacing = 0.75;
  int num_polar = 10; // 2

  double tolerance = 1e-5;
  int max_iters = 1000; //FIXME: 27 
  
  /* Create CMFD lattice */
  Cmfd* cmfd = new Cmfd();
  cmfd->useAxialInterpolation(false);

  std::vector<std::vector<double> > cmfd_widths{{0.05,0.63,0.63,0.63,0.63, 0.05},
                                                {0.05,0.63,0.63,0.63,0.63, 0.05},
                                                {1.25,1.25}};
  cmfd->setWidths(cmfd_widths);
  cmfd->setCentroidUpdateOn(true); 
  std::vector<std::vector<int> > cmfd_group_structure{{1,2,3},{4,5},{6,7}};
  cmfd->setGroupStructure(cmfd_group_structure);
  cmfd->setKNearest(1);


  /* Load the geometry */
  log_printf(NORMAL, "Creating geometry...");
  Geometry *geometry = new Geometry();
  geometry->loadFromFile(file, false); 
  
  geometry->setCmfd(cmfd);
#ifdef MPIx
  geometry->setDomainDecomposition(2, 2, 2, MPI_COMM_WORLD);
#endif
  geometry->initializeFlatSourceRegions();

  /* Create the track generator */
  log_printf(NORMAL, "Initializing the track generator...");
  Quadrature* quad = new EqualWeightPolarQuad();
  quad->setNumPolarAngles(num_polar);
  TrackGenerator3D track_generator(geometry, num_azim, num_polar, azim_spacing,
                                   polar_spacing);
  track_generator.setNumThreads(num_threads);
  track_generator.setQuadrature(quad);
  track_generator.setSegmentFormation(OTF_STACKS);
  std::vector<double> segmentation_zones {0.0, 1.0, 2.5};
  track_generator.setSegmentationZones(segmentation_zones);
  track_generator.generateTracks();

  /* Run simulation */
  CPULSSolver solver(&track_generator); 
  solver.setNumThreads(num_threads);
  solver.setVerboseIterationReport();
  solver.setConvergenceThreshold(tolerance);
  solver.computeEigenvalue(max_iters);
  solver.printTimerReport();
  double end_time = omp_get_wtime();
  log_printf(NORMAL, "Total Time = %6.4f", end_time - start_time);


  log_printf(TITLE, "Finished");
#ifdef MPIx
  MPI_Finalize();
#endif
  return 0;
}
