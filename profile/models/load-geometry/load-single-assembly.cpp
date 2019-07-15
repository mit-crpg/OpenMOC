#include "../../../src/CPUSolver.h"
#include "../../../src/CPULSSolver.h"
#include "../../../src/log.h"
#include "../../../src/Mesh.h"
#include <array>
#include <iostream>
#include "helper-code/group-structures.h"
#include <fenv.h>
#include <sstream>

int main(int argc,  char* argv[]) {

#ifdef MPIx
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED, &provided);
  log_set_ranks(MPI_COMM_WORLD);
#endif

  //feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT & ~FE_UNDERFLOW);
  /* Define geometry to load */
  //std::string file = "tw-single-assembly.geo";
  
  //std::string file = "tw-final-single-domain.geo";
  //std::string file = "tw-trunc-single-assembly.geo";
  std::string file = "tw-final-single-assembly.geo";

  /* FIXME
  int lat = 1;  
  std::string file = "lattice/tw-lattice-";
  std::stringstream ss;
  ss << lat;
  std::string lat_string = ss.str();
  file += lat_string;
  file += "x";
  file += lat_string;
  file += "x";
  file += lat_string;
  file += ".geo";
  */

  /* Define simulation parameters */
#ifdef OPENMP
  int num_threads = omp_get_num_threads();
#else
  int num_threads = 1;
#endif

  double azim_spacing = 0.05;
  int num_azim = 8; // 32
  double polar_spacing = 0.75;
  int num_polar = 4; // 10

  double tolerance = 1e-10;
  int max_iters = 50;

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
  double min_z = geometry.getMinZ();
  double max_z = geometry.getMaxZ();
  log_printf(NORMAL, "Height = %8.6e", max_z - min_z);
#ifdef MPIx
  //geometry.setDomainDecomposition(lat, lat, lat, MPI_COMM_WORLD); //FIXME 17x17xN
  geometry.setDomainDecomposition(1, 1, 100, MPI_COMM_WORLD); //FIXME 17x17xN
#endif
  //geometry.setNumDomainModules(1, 1, refines);
  geometry.setOverlaidMesh(2.0);
  geometry.initializeFlatSourceRegions();

  /* Create the track generator */
  log_printf(NORMAL, "Initializing the track generator...");
  TrackGenerator3D track_generator(&geometry, num_azim, num_polar, azim_spacing,
                                   polar_spacing);
  track_generator.setNumThreads(num_threads);
  //track_generator.setNumThreads(omp_get_max_threads());
  track_generator.setSegmentFormation(OTF_STACKS);
  double z_arr[] = {20., 34., 36., 38., 40., 98., 104., 150., 156., 202., 208.,
                    254., 260., 306., 312., 360., 366., 400., 402., 412., 414.,
                    418., 420.};
  //double z_arr[] = {120., 140.};
  //FIXME double z_arr[] = {min_z, max_z};
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
  solver.stabilizeTransport(0.25);
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
