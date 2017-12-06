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
  MPI_Init(&argc, &argv);
  log_set_ranks(MPI_COMM_WORLD);
#endif

  //feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT & ~FE_UNDERFLOW);
  /* Define geometry to load */
  //std::string file = "single-assembly-nc.geo";
  //std::string file = "single-assembly-5G-adjust.geo";
  //std::string file = "single-assembly-5G-adjust-no-neg.geo";
  //std::string file = "old-geo-files/single-no-TC.geo";
  //std::string file = "v1-single-assembly.geo";
  //std::string file = "v1-single-assembly-2D.geo";
  //std::string file = "final-single-assembly.geo";
  //std::string file = "assembly-lattice-1x1.geo";
  //std::string file = "tc-sa.geo";
  //std::string file = "rad-mesh/sa-fuel-4-mod-0.geo";
  //std::string file = "tw-single-assembly.geo";
  std::string file = "radial-mesh-geo/tw-single-assembly-"
                     "fr0fs4-mr1ms8-gr0gs8.geo";
  //std::string file = "tw-trunc-single-assembly.geo";

  /* Define simulation parameters */
  #ifdef OPENMP
  int num_threads = omp_get_max_threads();
  #else
  int num_threads = 1;
  #endif

  double azim_spacing = 0.05;
  int num_azim = 64;
  double polar_spacing = 0.75;
  int num_polar = 10;

  double tolerance = 1e-4;
  int max_iters = 500;

  /* Create CMFD lattice */
  Cmfd cmfd;
  cmfd.useAxialInterpolation(true);
  cmfd.setLatticeStructure(17, 17, 200); //FIXME 5 / 200
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
  geometry.loadFromFile(file, true);
  geometry.setCmfd(&cmfd); //FIXME OFF /ON
  log_printf(NORMAL, "Pitch = %8.6e", geometry.getMaxX() - geometry.getMinX());
  log_printf(NORMAL, "Height = %8.6e", geometry.getMaxZ() - geometry.getMinZ());
#ifdef MPIx
  geometry.setDomainDecomposition(1, 1, 100, MPI_COMM_WORLD); //FIXME 17x17xN
#endif
  //geometry.setNumDomainModules(17, 17, 1);
  geometry.setOverlaidMesh(2.0);
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
/* FIXME
  double z_arr[] = {120., 130.};
                    */
  //double z_arr[] = {0., 20., 380., 400.};
  std::vector<double> segmentation_zones(z_arr, z_arr + sizeof(z_arr)
                                         / sizeof(double));
  track_generator.setSegmentationZones(segmentation_zones);
  track_generator.generateTracks();

  /* Run simulation */
  CPULSSolver solver(&track_generator); //FIXME LS / FS
  solver.setNumThreads(num_threads);
  solver.setVerboseIterationReport();
  solver.setConvergenceThreshold(tolerance);
  solver.setCheckXSLogLevel(INFO);
  //solver.loadInitialFSRFluxes("fluxes-70-group");
  //solver.setResidualByReference("ref-fluxes-70-group-flat");
  //solver.setResidualByReference("ref-fluxes-70-group-flat-trunc");
  //solver.stabalizeTransport(1.0, YAMAMOTO);
  //solver.stabalizeTransport(1.0/16.0);
  solver.stabalizeTransport(0.25);
  solver.computeEigenvalue(max_iters);
  solver.printTimerReport();

  Lattice mesh_lattice;
  Mesh mesh(&solver);
  mesh.createLattice(17, 17, 200); //FIXME 1
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
