#include "../../../src/CPUSolver.h"
#include "../../../src/CPULSSolver.h"
#include "../../../src/log.h"
#include "../../../src/Mesh.h"
#include <array>
#include <iostream>
#include "helper-code/group-structures.h"
#include <fenv.h>

int main(int argc,  char* argv[]) {

#ifdef MPIx
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED, &provided);
  log_set_ranks(MPI_COMM_WORLD);
#endif

  /* Define geometry to load */
  //std::string file = "v11-full-core-new-xs.geo";
  //std::string file = "will-3D-TC-core.geo";
  //std::string file = "b2-full-core.geo";
  //std::string file = "final-full-core.geo";
  //std::string file = "final-full-core-full-TC.geo";
  //std::string file = "beavrs-casmo-water-xs.geo";
  std::string file = "tw-final-full-core-3D.geo";
  //feenableexcept(FE_DIVBYZERO | FE_OVERFLOW | FE_INVALID);
  //feenableexcept(FE_DIVBYZERO | FE_OVERFLOW | FE_INVALID);
  
  double start_time = omp_get_wtime();

  /* Define simulation parameters */
#ifdef OPENMP
  int num_threads = omp_get_num_threads();
#else
  int num_threads = 1;
#endif
 
  double azim_spacing = 0.05;
  int num_azim = 8; // 8
  double polar_spacing = 0.75;
  int num_polar = 2; // 2

  double tolerance = 1e-4;
  int max_iters = 27; //FIXME: 27
  
  /* Create CMFD lattice */
  Cmfd cmfd;
  cmfd.useAxialInterpolation(true);
  cmfd.setLatticeStructure(17*17, 17*17, 200);
  cmfd.setKNearest(1);
  std::vector<std::vector<int> > cmfd_group_structure =
      get_group_structure(70, 8);
  cmfd.setGroupStructure(cmfd_group_structure);
  cmfd.setCMFDRelaxationFactor(0.5);
  cmfd.setSORRelaxationFactor(1.6);
  cmfd.useFluxLimiting(true);
  //cmfd.rebalanceSigmaT(true);

  /* Load the geometry */
  log_printf(NORMAL, "Creating geometry...");
  Geometry geometry;
  geometry.loadFromFile(file, true);

  //FIXME
  int num_rad_discr = 96;
  int rad_discr_domains[96*2];
  int ind = 0;
  for (int i=0; i < 17; i++) {
    int offset = std::abs(8-i);
    if (offset == 8) {
      for (int j=0; j < 17; j++) {
        rad_discr_domains[2*ind] = i;
        rad_discr_domains[2*ind+1] = j;
        ind++;
      }
    }
    else if (offset == 7) {
      for (int j=0; j < 5; j++) {
        rad_discr_domains[2*ind] = i;
        rad_discr_domains[2*ind+1] = j;
        ind++;
      }
      for (int j=12; j < 17; j++) {
        rad_discr_domains[2*ind] = i;
        rad_discr_domains[2*ind+1] = j;
        ind++;
      }
    }
    else if (offset == 6) {
      for (int j=0; j < 3; j++) {
        rad_discr_domains[2*ind] = i;
        rad_discr_domains[2*ind+1] = j;
        ind++;
      }
      for (int j=14; j < 17; j++) {
        rad_discr_domains[2*ind] = i;
        rad_discr_domains[2*ind+1] = j;
        ind++;
      }
    }
    else if (offset == 5 || offset == 4) {
      for (int j=0; j < 2; j++) {
        rad_discr_domains[2*ind] = i;
        rad_discr_domains[2*ind+1] = j;
        ind++;
      }
      for (int j=15; j < 17; j++) {
        rad_discr_domains[2*ind] = i;
        rad_discr_domains[2*ind+1] = j;
        ind++;
      }
    }
    else if (offset < 4) {
      rad_discr_domains[2*ind] = i;
      rad_discr_domains[2*ind+1] = 0;
      ind++;
      rad_discr_domains[2*ind] = i;
      rad_discr_domains[2*ind+1] = 16;
      ind++;
    }
  }

  geometry.setCmfd(&cmfd); //FIXME OFF /ON
  log_printf(NORMAL, "Pitch = %8.6e", geometry.getMaxX() - geometry.getMinX());
  log_printf(NORMAL, "Height = %8.6e", geometry.getMaxZ() - geometry.getMinZ());
#ifdef MPIx
  //geometry.setDomainDecomposition(17, 17, 40, MPI_COMM_WORLD);
  //geometry.setNumDomainModules(1, 1, 1);
  //geometry.setDomainDecomposition(17, 17, 20, MPI_COMM_WORLD);
  //geometry.setNumDomainModules(1, 1, 2);
  geometry.setDomainDecomposition(17, 17, 5, MPI_COMM_WORLD);
  geometry.setNumDomainModules(1, 1, 8);
#endif
  geometry.setOverlaidMesh(2.0, 17*17*3, 17*17*3, num_rad_discr, 
                           rad_discr_domains);
  geometry.initializeFlatSourceRegions();

  /* Create the track generator */
  log_printf(NORMAL, "Initializing the track generator...");
  TrackGenerator3D track_generator(&geometry, num_azim, num_polar, azim_spacing,
                                   polar_spacing);
  track_generator.setNumThreads(num_threads);
  track_generator.setSegmentFormation(OTF_STACKS);
  double z_arr[] = {20., 34., 36., 38., 40., 98., 104., 150., 156., 202., 208.,
                    254., 260., 306., 312., 360., 366., 400., 402., 412., 414.,
                    418., 420.};
  std::vector<double> segmentation_zones(z_arr, z_arr + sizeof(z_arr) 
                                         / sizeof(double));
  track_generator.setSegmentationZones(segmentation_zones);
  track_generator.generateTracks();

  /* Run simulation */
  CPULSSolver solver(&track_generator);
  solver.stabilizeTransport(0.25);
  solver.setNumThreads(num_threads);
  solver.setVerboseIterationReport();
  solver.setConvergenceThreshold(tolerance);
  solver.setCheckXSLogLevel(INFO);
  solver.computeEigenvalue(max_iters);
  solver.printTimerReport();
  double end_time = omp_get_wtime();
  log_printf(NORMAL, "Total Time = %6.4f", end_time - start_time);

  Lattice mesh_lattice;
  Mesh mesh(&solver);
  mesh.createLattice(17*17, 17*17, 200);
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
